function m = measurable(hm, X, opts)
% FORWARD.MEASURABLE  What fraction of a HYPOTHESISED cortical pattern the array can see.
%
%   m = rheome.forward.measurable(hm, X)
%   m = rheome.forward.measurable(hm, X, opts)
%
% The whitened forward operator has a row space of dimension at most nCh (~270 for MEG) inside
% a source space of ~20k vertices. Everything orthogonal to it produces LITERALLY NO FIELD:
% add it to a source and the sensors do not move. This returns, for each column of X, the
% fraction of its energy that lies in that row space --
%
%     fraction = || V_r' * u ||^2 / || u ||^2 ,   u = sqrt(area) .* x,  V_r from svd(W*Gc*diag(sqrt(area)))
%
% -- so 0.9 means the data could in principle constrain 90% of that pattern and 0.1 means the
% pattern is nine-tenths invisible no matter how clean the recording.
%
% ⚠⚠ DO NOT RUN THIS ON YOUR OWN RECONSTRUCTION. IT RETURNS 1 BY CONSTRUCTION. A minimum-norm
% estimate is, by definition, the solution with the SMALLEST norm consistent with the data --
% which is the one whose null-space component is exactly ZERO. Every distributed linear inverse
% x = K*b likewise lives in range(K), which is the row space of the gain up to the prior's
% weighting. So MNE output scores ~1.0 always, and that number means nothing about whether the
% reconstruction is right. tMeasurable pins this as a test precisely so nobody reads it as a
% quality metric. What the null space actually means is the opposite: it is the set of sources
% you could ADD to the estimate without changing one sample of the data, so the reconstruction
% is one arbitrary representative of a huge equivalence class.
%
% ⭐ THE TWO USES THAT ARE NOT TAUTOLOGICAL, both on patterns you have hypothesised rather than
% reconstructed:
%
%   1. IS THIS HYPOTHESIS TESTABLE AT ALL? Build the pattern you are claiming -- a rotor at a
%      site, an ROI activation, an eigenmode -- and read its fraction. A low value says the
%      recording can neither confirm nor refute it, which is a result to report rather than a
%      reason to keep analysing.
%
%   2. ARE TWO HYPOTHESES DISTINGUISHABLE? What separates them is their DIFFERENCE, so pass
%      the difference: rheome.forward.measurable(hm, x1 - x2). Clockwise versus anticlockwise rotor at
%      the same site, or a rotor at site A versus site B -- if little of the difference is
%      measurable, no inverse and no statistic can tell them apart, because the two produce
%      near-identical sensor data.
%
% ⚠ RANK IS A CHOICE, AND THE ANSWER MOVES WITH IT. r is set by a relative singular-value cut,
% and rheome.forward.observability measured 94 / 183 / 255 modes above 1e-2 / 1e-3 / 1e-4 on a
% 270-channel MEG array -- so a pattern's "measurable fraction" is really "measurable at this
% assumed SNR". Report .rank alongside .fraction, and sweep it if a conclusion depends on it.
%
% INPUTS:
%   hm    rheome.io.read.headmodel struct (.Gain [nCh x 3nV], .GridOrient [nV x 3], .nV)
%   X     [nV x nCol]  scalar current DENSITY along the constrained orientation, or
%         [3nV x nCol] unconstrained current density (xyz interleaved per vertex)
%   opts  .Channels    gain rows to keep (default all)
%         .NoiseCov    [nCh x nCh] over those rows, [] for identity
%         .Reg         loading as a fraction of mean(diag(C))   (default 0.05)
%         .SVCut       relative singular-value cut setting r    (default 1e-3)
%         .Area        [nV x 1] vertex areas (default: ones -- pass sum(Mass,2) to weight properly)
%         .Project     true to also return the measurable part of X   (default false)
%
% OUTPUT (struct m):
%   .fraction [nCol x 1]   measurable energy fraction, in [0,1]
%   .rank                  r actually used;  .sv  singular values;  .SVCut
%   .projected [same size as X]  the measurable part, if opts.Project
%
% See also: rheome.forward.observability, rheome.flow.crosstalk, rheome.flow.sensitivity
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(opts), opts = struct(); end
    ch    = i_opt(opts, 'Channels', []);
    C     = i_opt(opts, 'NoiseCov', []);
    reg   = i_opt(opts, 'Reg',      0.05);
    cut   = i_opt(opts, 'SVCut',    1e-3);
    doPrj = i_opt(opts, 'Project',  false);

    G  = double(hm.Gain);  nV = hm.nV;
    if isempty(ch), ch = 1:size(G,1); end
    G  = G(ch(:).', :);  nCh = size(G,1);

    area = i_opt(opts, 'Area', ones(nV,1));
    area = double(area(:));
    if numel(area) ~= nV
        error('forward:measurable:area', 'Area has %d entries, expected %d.', numel(area), nV);
    end
    s = sqrt(max(area, realmin));

    if size(X,1) == nV
        O = double(hm.GridOrient);
        if isempty(O)
            error('forward:measurable:orient', 'Scalar X needs hm.GridOrient.');
        end
        O = O ./ max(vecnorm(O,2,2), eps);
        Osp = sparse((1:3*nV).', repelem((1:nV).',3,1), reshape(O.',[],1), 3*nV, nV);
        Go  = G * Osp;                         % [nCh x nV]
        sw  = s;
    elseif size(X,1) == 3*nV
        Go  = G;                               % [nCh x 3nV]
        sw  = repelem(s, 3, 1);
    else
        error('forward:measurable:size', ...
            'X has %d rows; expected %d (constrained) or %d (unconstrained).', size(X,1), nV, 3*nV);
    end

    % --- whitener over the SAME rows
    if isempty(C)
        W = speye(nCh);
    else
        C = double(C);
        if ~isequal(size(C), [nCh nCh])
            error('forward:measurable:cov', 'NoiseCov is %s but %d channels selected.', ...
                mat2str(size(C)), nCh);
        end
        C = (C + C.')/2 + reg * mean(diag(C)) * eye(nCh);
        [U,d] = eig(C, 'vector');
        W = diag(1./sqrt(max(d, max(d)*1e-12))) * U.';
    end

    % --- row space in mass-orthonormal source coordinates
    Ao = W * (Go .* sw.');                     % [nCh x n]  (scales columns by sqrt(area))
    [~, sv, Vs] = svd(Ao, 'econ');
    sv = diag(sv);
    r  = max(1, sum(sv > sv(1)*cut));
    Vr = Vs(:, 1:r);

    Uu = sw .* double(X);                      % source coords with unit mass metric
    Pu = Vr.' * Uu;                            % [r x nCol]
    den = sum(Uu.^2, 1).';
    m.fraction = sum(Pu.^2, 1).' ./ max(den, realmin);
    m.fraction(den <= 0) = NaN;                % an all-zero column has no fraction to report
    m.rank = r;  m.sv = sv;  m.SVCut = cut;  m.nCh = nCh;

    if doPrj
        m.projected = (Vr * Pu) ./ sw;
    end
end

function v = i_opt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
