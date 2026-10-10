function ctx = context(name, K, Klbo, opts)
% FLOW.CONTEXT  One-time precompute bundle for the fused flow kernels.
%   ctx = rheome.flow.context(name [,K [,Klbo [,opts]]])
%   ctx = rheome.flow.context(name, [], Klbo, struct('Method','dirac'))
%
% ⭐⭐ THE SOURCE ESTIMATE IS A PLAIN MINIMUM NORM BY DEFAULT (opts.Method = 'mne'). The flow
% quantities are differential operators applied to a current field, and the most basic current
% field is the raw whitened minimum-norm estimate on the full unconstrained leadfield. The LBO
% eigenbasis then enters AFTERWARDS, on the resulting scalar maps, purely as the analysis tool
% that separates scales and sets the observation aperture. Those are two separate jobs and this
% function keeps them separate.
%
% ⚠ THIS DEFAULT CHANGED, AND EVERY DOWNSTREAM NUMBER MOVES WITH IT. It used to be the Dirac
% route unconditionally. That put a SECOND, INVISIBLE band limit ahead of the analysis one:
% projecting the leadfield onto K Dirac modes truncates the current before the inverse is even
% solved. Measured on a reference subject at the pipeline defaults (400 modes/hemisphere both sides):
%
%     LBO analysis axis reaches   28.5 mm   (frame floor sigma >= 14.0 mm, printed below)
%     Dirac current stops at      90.4 mm
%     -> 3.2x mismatch
%
% ⚠ "STOPS AT" DEPENDS ON THE DEFINITION, and the two in use here differ by 1.7x. The 90.4 mm
% above was recorded without one; a per-mode calibration of the SAME cached basis
% (project each Dirac mode's three ambient current components onto
% the LBO basis, take the energy-weighted wavenumber centroid) puts the FINEST mode's own scale
% at 52 mm, with 396 of the 800 modes resolving to one hemisphere and spanning 459 down to
% 52 mm. Those are different quantities -- the finest mode's scale is not the scale at which a
% RECONSTRUCTION from all of them still carries fidelity, which is coarser and is plausibly
% what 90.4 mm measured. Quote whichever the claim needs, and say which it is: a mismatch ratio
% computed against the wrong one is wrong by that 1.7x.
%
% ⚠ AND THE DIRAC SPECTRUM CARRIES NO LENGTH AXIS OF ITS OWN. Its eigenvalues run 0 to 9.5e-6
% in the relative operator's units, arrive in groups of four (the quaternion structure), and
% 2*pi/sqrt(lambda) returns millions of millimetres; dbasis.Scales is [2 x 2], per-hemisphere
% normalisation, not per mode. graphfilterbank therefore cannot be pointed at this basis --
% it assumes the operator's own lambda IS the scale axis. Calibrate first.
%
% so every aperture finer than ~90 mm was filtering structure the source estimate could not
% represent, and nothing in the pipeline said so. 'mne' has no such hidden floor: the source
% space is the full 3nV and the only band limit left is the analysis one, which IS reported.
% Pass opts.Method = 'dirac' for the old behaviour -- it is a legitimate alternative, with a
% geometry-aware prior, but it is not the baseline and it must be asked for.
%
% ⚠ K IS THE DIRAC MODE COUNT (Method='dirac' only); Klbo IS THE LBO ONE. They are independent,
% and it is Klbo that sets
% the lambda axis every scale, size and speed is read off. The method's default is Klbo = 400
% modes per hemisphere (Ks = 800), the value the pipeline and the paper use: pass it. Leaving Klbo
% empty takes whatever the richest cached LBO basis happens to be -- fine while exploring, NOT
% fine for a result, because the thing that changes is sigma and nothing announces it. K = 400
% when empty is the DIRAC mode count, not this one.
%
% The field kernel currentKernel [3V x C] is the source estimate in AMPLITUDE measure (physical
% current A.m, the right quantity for flow): for Method = 'mne' (default) the whitened minimum norm
% of rheome.inverse.mne on the full unconstrained leadfield; for Method = 'dirac' the Dirac source
% mapping of rheome.source.dirac, reconstructed as rheome.forward.reconstruct(diracInverse, dbasis).
% It then assembles the shared operators (face_gradient) and the Laplace–Beltrami eigenbasis used
% by the potentials. The LBO eigenbasis is LOADED from the cached per-hemisphere bases
% (rheome.load.bases) and assembled block-diagonally -- NO eigensolve. `name` is a cached dataset
% name (see rheome.import.study, rheome.import.bases; rheome.import.dataset for 'dirac'). Build
% this ONCE and reuse it across the rheome.flow.* builders.
%
% OUTPUT (struct ctx):
%   .currentKernel [3V x C]  the source estimate every flow kernel is built from (amplitude, A.m)
%   .Method  'mne' | 'dirac'                         .ImagingKernel [3V x C] same, as returned
%   .diracInverse [P x C] / .dbasis                  EMPTY unless Method = 'dirac'
%   .S surface   .fg face_gradient bundle            .lbo (.Phi [V x Ks], .Lambda [Ks x 1])
%   .F [C x nT] selected MEG data   .sfreq   .iSel selected channel rows   .chNames {1 x C}
%   .V nVert   .C nChan(sel)   .P nModes   .Ks LBO modes used by the potentials
%
% See also: rheome.inverse.mne, rheome.source.dirac, rheome.inverse.dirac, rheome.forward.reconstruct, rheome.operators.face_gradient
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(K),   K    = 400;       end
    if nargin < 3,                 Klbo = [];        end
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'Method') || isempty(opts.Method), opts.Method = 'mne'; end
    method = lower(char(opts.Method));

    ctx = struct();
    ctx.Method = method;
    switch method
        case 'mne'
            % ⭐ NO rheome.source.dirac AND NO DIRAC BASIS ON THIS PATH. Besides being the point, it
            % is why this branch is much lighter: the 369 MB Dirac basis is never read.
            st = rheome.load.study(name);
            S  = rheome.load.surface(name);
            if st.hm.nV ~= S.nV
                error('flow:context:grid', ...
                    'Leadfield grid has %d sources but the cortex has %d vertices.', st.hm.nV, S.nV);
            end
            isMEG = strcmpi(st.chan.Type, 'MEG');
            iSel  = find(isMEG(:) & (st.rec.ChannelFlag(:) == 1));   % same selection as rheome.source.dirac
            % ⚠ TAKE WHAT IS NEEDED AND DROP THE STUDY BEFORE THE SVD. st holds the whole
            % recording (~1 GB on a resting run); keeping it alive across an [nCh x 3nV]
            % decomposition is what turns a comfortable build into a swapping one.
            Gsel  = st.hm.Gain(iSel,:);
            ncm   = struct('NoiseCov',     st.ncov.NoiseCov(iSel,iSel), ...
                           'FourthMoment', i_sub(st.ncov.FourthMoment, iSel), ...
                           'nSamples',     i_sub(st.ncov.nSamples, iSel));
            chTypes = st.chan.Type(iSel);
            F = []; if ~isempty(st.rec.F), F = st.rec.F(iSel,:); end
            sfreq = st.rec.sfreq;
            clear st
            Res = rheome.inverse.mne(Gsel, ncm, ...
                    struct('ChannelTypes', {chTypes}, 'InverseMeasure', 'amplitude', ...
                           'nVert', S.nV));   % closes the divisible-by-3 ambiguity
            clear Gsel
            ctx.diracInverse  = [];
            ctx.dbasis        = [];
            ctx.ImagingKernel = Res.ImagingKernel;
            ctx.currentKernel = Res.ImagingKernel;                  % [3V x C]
            ctx.inverse       = Res;
            out = struct('S', S, 'F', F, 'sfreq', sfreq, 'iSel', iSel);
        case 'dirac'
            out = rheome.source.dirac(name, [], [], K, 'amplitude');   % amplitude = physical current (A.m)
            Res = out.Results;
            ctx.diracInverse  = Res.ImagingKernelMode;                      % [P x C]
            ctx.ImagingKernel = Res.ImagingKernel;                          % [3V x C] (amplitude)
            ctx.dbasis        = out.dbasis;
            ctx.currentKernel = rheome.forward.reconstruct(ctx.diracInverse, ctx.dbasis);    % [3V x C]
            ctx.inverse       = Res;
        otherwise
            error('flow:context:method', ...
                'Method must be ''mne'' (default) or ''dirac'', got ''%s''.', opts.Method);
    end
    ctx.S             = out.S;
    ctx.fg            = rheome.operators.face_gradient(out.S.Vertices, out.S.Faces);

    % Laplace-Beltrami eigenbasis for the flow coefficients + Helmholtz potentials. LOAD the cached
    % PER-HEMISPHERE bases (rheome.import.bases: K modes/hemisphere -- the optimal resolution -- so 2K whole-
    % brain) and assemble the whole-cortex basis block-diagonally. NO eigensolve on the analysis path:
    % a whole-brain solve would be both redundant AND ill-conditioned, since the disconnected cortex
    % has a 2-D null space (one constant per hemisphere) that the shift-invert solver chokes on.
    [B, binfo] = rheome.load.bases(name, 0.5, K, Klbo);
    labs = intersect(fieldnames(B), {'L','R'}, 'stable');
    Ktot = sum(cellfun(@(l) numel(B.(l).lbo.Lambda), labs));
    Phi  = zeros(ctx.S.nV, Ktot);  Lambda = zeros(Ktot, 1);  Mass = sparse(ctx.S.nV, ctx.S.nV);  col = 0;
    for i = 1:numel(labs)
        l = labs{i};  gv = double(B.(l).gv(:));  kh = numel(B.(l).lbo.Lambda);
        Phi(gv, col+(1:kh)) = B.(l).lbo.Phi;
        Lambda(col+(1:kh))  = B.(l).lbo.Lambda(:);
        Mass(gv, gv)        = B.(l).lbo.Mass;
        col = col + kh;
    end
    ctx.lbo = struct('Phi', Phi, 'Lambda', Lambda, 'Mass', Mass);
    ctx.Ks  = Ktot;
    % Say out loud which lambda axis this context is on, and what it costs at the fine end. The
    % frame's finest usable member is sigma_m >= 3.08/sqrt(lambda_max); a detection below that is a
    % floor rather than a measurement, and it is invisible unless it is printed.
    ctx.bases = binfo;
    fprintf('rheome.flow.context[%s]: inverse = %s | LBO K = %d/hemisphere (%d total) | lambda_max = %.4g | frame floor sigma_m >= %.1f mm\n', ...
        char(name), upper(method), binfo.K, Ktot, binfo.lambdaMax, 1000*binfo.sigmaFloor);

    ctx.F             = out.F;                                       % [C x nT] already selected
    ctx.sfreq         = out.sfreq;
    ctx.iSel          = out.iSel(:)';
    ctx.chNames       = arrayfun(@(i) sprintf('MEG%d', i), ctx.iSel, 'UniformOutput', false);
    ctx.V             = ctx.fg.nV;
    ctx.C             = size(ctx.currentKernel, 2);
    % ⚠ P IS THE SOURCE-SPACE DIMENSION THE INVERSE SOLVED IN, not a mode count. For 'mne' that
    % is the full 3nV -- there is no source-side truncation to report, which is the point.
    ctx.P             = size(ctx.currentKernel, 1);
    if strcmp(method,'dirac'), ctx.P = size(ctx.diracInverse, 1); end
end

function A = i_sub(A, idx)
    if ~isempty(A) && size(A,1) >= max(idx), A = A(idx, idx); end
end

% Author: Diellor Basha, 2026
