function out = ridges(f, S, basis, sigmas, opts)
% DETECT.RIDGES  Scale-space ridge tracking: detection, size and continuity from one volume.
%
%   out = rheome.detect.ridges(f, S, basis, sigmas)        % legacy: hand-declared sigma list (m)
%   out = rheome.detect.ridges(f, S, basis, frame)         % preferred: a designed frame (rheome.filters.frame)
%   out = rheome.detect.ridges(f, S, basis, ..., opts)
%
% Builds the scale-space volume  W(v,m,t) = |mexhat_{sigma_m}(Delta) f(:,t)|,  takes local maxima
% JOINTLY in (vertex, scale) -- so each detection carries a POSITION and a SIZE from one operation --
% and links them through time by ridge connectivity. The linking gate is the MESH NEIGHBOURHOOD
% itself: there is no tuned association distance, no lost-track buffer and no merge/split bookkeeping.
%
% WHY. Detect-then-associate discards continuity when each frame is reduced to a point list, then
% reconstructs it by matching. A ridge never discards it. Measured consequence on real alpha:
% rheome.detect.track caps motif lifetime at the minimum-track-length floor regardless of band, whereas
% ridges on the demodulated envelope reach ~1.3 s (see the grant workspace notes).
%
% ⚠ RUN THIS ON A PHASE-INVARIANT FIELD. An oscillating curl/divergence field reverses sign every
% half cycle, driving |W| to zero twice per period, so ridges terminate at each reversal exactly as
% tracks do. Pass the analytic ENVELOPE (demodulated), not the raw band-limited field, unless you
% specifically want the within-half-cycle structure.
%
% NORMALISATION. rheome.filters.mexhat(Lambda, t) = (t*lambda).*exp(-t*lambda) IS Lindeberg's
% scale-normalised Laplacian t*grad^2 G, so its PEAK response is scale-independent (e^-1 for every t)
% and the scale axis is directly comparable -- which is what makes the joint (vertex,scale) maximum a
% valid size estimate.
%
% ⚠ .sigma AND .sigmaMode ARE MEMBER LABELS, NOT VORTEX SIZES -- multiply by sqrt(2). The member
% whose gain PEAKS with a sigma-vortex's vorticity spectrum is the one with sigma_m = sigma, but the
% detector maximises the RESPONSE, which integrates the overlap against the mode density (uniform in
% lambda on a surface, by Weyl). That integral,
%
%     Resp(t) ~ int_0^inf (t*lam*e^{-t*lam})(lam*e^{-lam*sigma^2/2}) dlam = 2t/(t + sigma^2/2)^3
%
% is maximised at t = sigma^2/4, so the reported member sits at sigma/sqrt(2):
%
%     ⭐ gamma = sigma_m / sigma_true = 1/sqrt(2) = 0.7071   ->   sigma_true = sqrt(2) * sigma_m
%
% Analytic, not fitted; confirmed at 0.735 over sigma = 8-40 mm in the sphere validation. It is
% NOT applied here, because doing so silently would change every already-reported size -- apply it
% at the point of reporting, and say that you have.
%
% ⚠ AND THE FINEST MEMBERS UNDER-RESPOND. rheome.filters.frame sets t_min = 1/lambda_max, which puts 73.6%
% of the finest member's mass beyond the end of the basis (and 45.9% of the next one's), so the
% scale-space maximum is pushed one member coarser for small vortices. Members with
% sigma_m >= 3.08/sqrt(lambda_max) lose under 5%; below that, treat the size as a floor.
%
% INPUTS:
%   f       [nV x nT] scalar field over time (real or complex; |.| is taken)
%   S       surface struct (.Vertices, .Faces, .nV)
%   basis   Laplace-Beltrami eigenbasis (.Phi, .Lambda, .Mass, .nV)
%   sigmas  [1 x M] spatial scales in METRES. Scales with 2/sigma^2 > max(Lambda) are REFUSED
%           (the mexhat peak would fall above the basis's spectral limit and the filter would
%           silently clamp to a spectral edge -- a real artifact source).
%   opts    .thresh    'mad' (default, 3*robust sd of |W|) | a numeric absolute threshold
%           .minFrames minimum ridge length in frames (default 3)
%           .nbr       precomputed neighbour cell array (saves time across repeated calls)
%
% OUTPUT (struct out):
%   .ridges  struct array: .v (vertex path) .m (scale index path) .t (frames)
%                          .pos [n x 3] .sigma [n x 1] (m) .w [n x 1] (response)
%                          .nFrames .lifetime .netDisplacement .pathLength .speedNet
%                          .sigmaMode (modal scale, m) .nScales (distinct scales visited)
%   .summary .nRidges .nDetPerFrame .sigmas .scaleHist .lifetime .speedNet .thresh .dt
%   .W       the scale-space volume [nV x M x nT], single (only if opts.keepVolume)
%
% See also: rheome.filters.mexhat, rheome.filters.apply, rheome.detect.track, rheome.show.ridges
%
% Author: Diellor Basha, 2026

    if nargin < 5, opts = struct; end
    if ~isfield(opts,'thresh')     || isempty(opts.thresh),     opts.thresh = 'mad'; end
    if ~isfield(opts,'threshFrac') || isempty(opts.threshFrac), opts.threshFrac = 0.15; end
    if ~isfield(opts,'minFrames')  || isempty(opts.minFrames),  opts.minFrames = 3;  end
    if ~isfield(opts,'dt')         || isempty(opts.dt),         opts.dt = 1;         end
    if ~isfield(opts,'keepVolume') || isempty(opts.keepVolume), opts.keepVolume = false; end

    nV = S.nV;  nT = size(f,2);
    lam = basis.Lambda(:);

    % ---- the scale bank: a designed FRAME, or a hand-declared sigma list (legacy) ----
    if isstruct(sigmas) && isfield(sigmas, 'g')
        % DESIGNED FRAME (rheome.filters.frame). Preferred: members are log-spaced from t_min = 1/lambda_max,
        % so every one is representable by construction and nothing has to be refused; a low-pass
        % scaling function catches energy coarser than the coarsest wavelet, which a bare sigma list
        % discards silently; and coverage is quantifiable via rheome.filters.frame_bounds.
        frame  = sigmas;
        G      = rheome.filters.frame_gains(frame, lam);                  % [Ks x M]
        M      = size(G,2);
        % Label members with the EXACT scale sqrt(2*t), not the gain-weighted centroid: the centroid
        % is contaminated by where the eigenvalue axis was truncated and compresses the fine end
        % badly (on a reference subject it labelled the two finest members 9.1 and 9.5 mm when they are 7.1 and
        % 9.6 mm -- two distinct scales reported as one).
        if isfield(frame,'Sigma') && ~isempty(frame.Sigma)
            sigmas = frame.Sigma;
        else
            sigmas = sqrt(2) ./ max(frame.Centers, eps);
        end
        if strcmpi(frame.Family,'mexhat'), sigmas(1) = Inf; end    % member 1 is the scaling function
        if isfield(frame,'Gamma'),  gamma  = frame.Gamma;  else, gamma  = 1/sqrt(2); end
        if isfield(frame,'Usable'), usable = frame.Usable; else, usable = true(1,M); end
        bnd = rheome.filters.frame_bounds(frame, lam);
        if bnd.A <= 0
            warning('detect:ridges:coverage', ...
                'Frame lower bound A = %.3g: part of the spectrum is uncovered (widen lrange).', bnd.A);
        end
    else
        % LEGACY sigma list: representability is not guaranteed, so unrepresentable scales are
        % REFUSED rather than silently clamped to something the basis cannot support.
        sigmas = sigmas(:)';
        rep = (2 ./ sigmas.^2) < max(lam);
        if ~all(rep)
            warning('detect:ridges:scale', ...
                'Refusing %d scale(s) with 2/sigma^2 > lambda_max = %.4g: %s mm', ...
                sum(~rep), max(lam), mat2str(1000*sigmas(~rep)));
            sigmas = sigmas(rep);
        end
        M = numel(sigmas);
        if M == 0, error('detect:ridges:noScales', 'No representable scales.'); end
        G = rheome.filters.mexhat(lam, sigmas.^2/2);        % [Ks x M], Lindeberg-normalised
        % Same members, same calibration. Representability was checked above (2/sigma^2 < lambda_max),
        % but that is a much weaker bar than "responds properly": a member whose peak sits just under
        % lambda_max still loses most of its mass past the end of the basis.
        gamma  = 1/sqrt(2);
        u      = (sigmas.^2/2) * max(lam);
        usable = (1 + u) .* exp(-u) <= 0.05;
        if ~all(usable)
            warning('detect:ridges:truncated', ...
                ['%d of %d scales lose >5%% of their mass beyond lambda_max: %s mm. They ' ...
                 'under-respond, so sizes near the fine end read too coarse.'], ...
                sum(~usable), M, mat2str(round(1000*sigmas(~usable),1)));
        end
    end

    % ---- mesh adjacency ----
    if isfield(opts,'nbr') && ~isempty(opts.nbr)
        nbr = opts.nbr;
    else
        F = S.Faces;
        A = sparse([F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1)], ...
                   [F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3)], 1, nV, nV) > 0;
        nbr = cell(nV,1);
        for v = 1:nV, nbr{v} = find(A(v,:)); end
    end

    % ---- 1-2. coefficients (shared) -> scale-space volume ----
    c = basis.Phi' * (basis.Mass * f);
    W = zeros(nV, M, nT, 'single');
    for m = 1:M
        W(:,m,:) = single(abs(basis.Phi * (G(:,m) .* c)));
    end

    % ---- threshold. A pure MAD rule fails on sparse fields: if most of the volume is ~0 (one
    %      synthetic structure on an otherwise empty sphere) then MAD -> 0 and everything above
    %      numerical noise becomes a detection. Floor it at a fraction of the volume's own peak, so
    %      the rule behaves the same on a clean synthetic field and on broadband real data. ----
    if ischar(opts.thresh) || isstring(opts.thresh)
        mad3 = 3 * 1.4826 * median(abs(W(:) - median(W(:))));
        thr  = max(mad3, opts.threshFrac * max(W(:)));
    else
        thr = opts.thresh;
    end

    % ---- 3. local maxima jointly in (vertex, scale) ----
    dV = cell(nT,1); dM = cell(nT,1); dW = cell(nT,1);
    for t = 1:nT
        Wt = W(:,:,t);
        vv = zeros(1,0); mm = zeros(1,0); ww = zeros(1,0);
        for m = 1:M
            ms = max(1,m-1):min(M,m+1);
            cand = find(Wt(:,m) > thr).';
            for v = cand
                loc = Wt([v nbr{v}], ms);
                if Wt(v,m) >= max(loc(:))
                    vv(end+1)=v; mm(end+1)=m; ww(end+1)=Wt(v,m); %#ok<AGROW>
                end
            end
        end
        dV{t}=vv; dM{t}=mm; dW{t}=ww;
    end

    % ---- 4. ridge linking: mesh neighbourhood + one scale step. No tuned gate. ----
    R = {};  active = zeros(1,0);
    for t = 1:nT
        v = dV{t}; m = dM{t}; w = dW{t};
        used = false(1,numel(v));  newActive = zeros(1,0);
        for a = active
            lv = R{a}.v(end); lm = R{a}.m(end);
            cand = find(~used & abs(m-lm) <= 1 & ismember(v, [lv nbr{lv}]));
            if isempty(cand), continue; end
            [~,b] = max(w(cand));  pick = cand(b);
            R{a}.v(end+1)=v(pick); R{a}.m(end+1)=m(pick); R{a}.t(end+1)=t; R{a}.w(end+1)=w(pick);
            used(pick) = true;  newActive(end+1) = a; %#ok<AGROW>
        end
        for k = find(~used)
            R{end+1} = struct('v',v(k),'m',m(k),'t',t,'w',w(k)); %#ok<AGROW>
            newActive(end+1) = numel(R); %#ok<AGROW>
        end
        active = newActive;
    end

    % ---- descriptors ----
    keep = cellfun(@(r) numel(r.t) >= opts.minFrames, R);
    R = R(keep);
    ridge = i_empty();
    for k = 1:numel(R)
        r = R{k};
        P = S.Vertices(r.v(:), :);
        d.v = r.v(:);  d.m = r.m(:);  d.t = r.t(:);  d.w = r.w(:);
        d.pos = P;  d.sigma = sigmas(r.m(:)).';
        d.nFrames = numel(r.t);
        d.lifetime = (r.t(end) - r.t(1)) * opts.dt;
        d.netDisplacement = norm(P(end,:) - P(1,:));
        d.pathLength = sum(vecnorm(diff(P,1,1), 2, 2));
        d.speedNet = d.netDisplacement / max(d.lifetime, eps);
        d.sigmaMode = sigmas(mode(r.m(:)));
        d.nScales = numel(unique(r.m));
        % The MEMBER label divided by the calibration gamma = 1/sqrt(2): the vortex size this
        % detection corresponds to. Carried alongside .sigma rather than replacing it, so a reader
        % can see which one a figure or table was built from.
        d.sigmaTrue     = d.sigma / gamma;
        d.sigmaModeTrue = d.sigmaMode / gamma;
        ridge(end+1) = d; %#ok<AGROW>
    end

    out.ridges = ridge;
    out.summary.nRidges      = numel(ridge);
    out.summary.nDetPerFrame = cellfun(@numel, dV);
    out.summary.sigmas       = sigmas;
    out.summary.gamma        = gamma;              % sigma_m / sigma_true = 1/sqrt(2), analytic
    out.summary.sigmasTrue   = sigmas / gamma;     % what the members correspond to as vortex sizes
    out.summary.usable       = usable;             % false where the member loses >5% past lambda_max

    % ⚠ IS THE BANK EVEN COVERING WHERE THE ENERGY IS? The mexhat scaling function is an unbounded
    % low-pass: it catches everything coarser than the coarsest wavelet, and its "size" is Inf. If
    % most detections land there, the motifs are COARSER than the bank reaches and no size has been
    % measured at all -- sigmaMode comes back Inf, which is not a number anyone should report.
    % This happens on raising K at fixed Nf: rheome.filters.frame ties t_max to lambda_max, so a richer
    % basis shifts the WHOLE bank finer and abandons the coarse end (measured on a reference subject: the
    % coarsest member went 46.6 mm at K=400 to 27.9 mm at K=1000, and the low-pass went from 1.7%
    % to 11.5% of the energy). The fix is more members, not more modes.
    if ~isempty(ridge)
        nLP = sum(~isfinite([ridge.sigmaMode]));
        out.summary.fracLowpass = nLP / numel(ridge);
        if out.summary.fracLowpass > 0.25
            warning('detect:ridges:lowpass', ...
                ['%.0f%% of ridges peak on the SCALING FUNCTION (an unbounded low-pass), so their ' ...
                 'size is Inf and unmeasured: the motifs are coarser than the coarsest wavelet ' ...
                 '(%.1f mm). Add frame members at the coarse end -- raising K alone makes this ' ...
                 'WORSE, since it shifts the whole bank finer.'], ...
                100*out.summary.fracLowpass, 1000*max(sigmas(isfinite(sigmas))));
        end
    else
        out.summary.fracLowpass = NaN;
    end
    out.summary.scaleHist    = histcounts([dM{:}], 0.5:1:(M+0.5));
    out.summary.lifetime     = [ridge.lifetime];
    out.summary.speedNet     = [ridge.speedNet];
    out.summary.thresh       = thr;
    out.summary.dt           = opts.dt;
    if opts.keepVolume, out.W = W; end
end

function d = i_empty()
    d = struct('v',{},'m',{},'t',{},'w',{},'pos',{},'sigma',{},'nFrames',{},'lifetime',{}, ...
               'netDisplacement',{},'pathLength',{},'speedNet',{},'sigmaMode',{},'nScales',{}, ...
               'sigmaTrue',{},'sigmaModeTrue',{});
end

% Author: Diellor Basha, 2026
