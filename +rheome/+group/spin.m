function R = spin(M, Ref, U, hemi, opts)
% GROUP.SPIN  Does an effect map follow reference maps more than its own spatial smoothness explains? Spin test, FWER over the references.
%
%   R = rheome.group.spin(M, Ref, U, hemi)
%   R = rheome.group.spin(M, Ref, U, hemi, NSpin=1000, Seed=1)
%   R = rheome.group.spin(M, Ref, U, hemi, Index=R0.index)       % the same rotations as an earlier call
%   R = rheome.group.spin(M, Ref, U, hemi, Cortex=G.cortex_area > 0)   % medial-wall tiles out
%
%   M     [m x 1] the tile map tested (e.g. rheome.group.glm's t)
%   Ref   [m x K] K reference maps on the same tiles
%   U     [m x 3] tile centres on the template sphere (unit vectors)
%   hemi  [m x 1] "L" / "R"
%
% r_k = Pearson(M, Ref_k). The null rotates M on the sphere (Alexander-Bloch et al. 2018): one random rotation
% per spin for both hemispheres, the right one mirrored (x -> -x), each tile taking the value of the tile nearest
% its rotated centre within its hemisphere, as rheome.modes.fit's N2. p_k = (1 + #{|r*_k| >= |r_k|}) / (1 + NSpin);
% p_fwer_k compares |r_k| with max_k |r*_k| of each spin (FWER over the K references).
%
% A NaN r (M or a reference constant, or holding NaN) gets p = p_fwer = 1.
%
% Returns R with .r .p .p_fwer [1 x K], .rnull [NSpin x K] and .index [NSpin x m] (spin s maps tile i to tile
% index(s, i); pass it back as Index= to test many maps against one set of rotations).
%
% ⚠ FEW TILES, FEW DISTINCT SPINS. Nearest-tile reassignment on m/2 tiles per hemisphere repeats maps across
% spins (m = 16 at depth 3) and many rotations land a tile on itself: ties are counted against the effect (>=),
% so the test is conservative there and cannot reach small p. Read its calibration (rheome.group.calibrate).
%
% MEDIAL WALL (decision 3882c4ce, 2026-10-14). Cortex= marks the tiles that hold cortex (rheome.geom.tiles
% .cortex_area > 0). The others are dropped before anything else: they are neither scored nor a target of a
% rotated tile (the rotated wall is discarded and reassigned to the nearest cortex tile), and .index then
% numbers the kept tiles (R.cortex says which).
%
% See also: rheome.group.glm, rheome.modes.fit
%
% Author: Diellor Basha, 2026

    arguments
        M (:,1) double
        Ref double
        U (:,3) double
        hemi (:,1) string
        opts.NSpin (1,1) double {mustBeInteger, mustBePositive} = 1000
        opts.Seed (1,1) double = 1
        opts.Index double = []
        opts.Cortex logical = logical.empty
    end
    if ~isempty(opts.Cortex) && ~all(opts.Cortex)
        k = opts.Cortex(:);
        R = rheome.group.spin(M(k), Ref(k, :), U(k, :), hemi(k), NSpin=opts.NSpin, Seed=opts.Seed, Index=opts.Index);
        R.cortex = k;  return
    end
    m = numel(M);
    if size(Ref, 1) ~= m || size(U, 1) ~= m || numel(hemi) ~= m, error('group:spin:size', 'M, Ref, U and hemi need one row per tile'); end
    R.index = opts.Index;
    if isempty(R.index), R.index = i_index(U, hemi, opts.NSpin, opts.Seed); end
    R.r = corr(M, Ref);  ar = abs(R.r);
    R.rnull = abs(corr(M(R.index'), Ref));                  % [NSpin x K]
    R.rnull(isnan(R.rnull)) = 0;                           % a spin that maps every tile to one value: no correlation
    ns = size(R.index, 1);
    R.p = (1 + sum(R.rnull >= ar - 1e-12, 1)) / (1 + ns);
    R.p_fwer = (1 + sum(max(R.rnull, [], 2) >= ar - 1e-12, 1)) / (1 + ns);
    R.p(isnan(R.r)) = 1;  R.p_fwer(isnan(R.r)) = 1;        % a NaN or constant map (or reference) is not evidence
end

function I = i_index(U, hemi, n, seed)
    U = U ./ vecnorm(U, 2, 2);  rs = RandStream('mt19937ar', 'Seed', seed);  Fx = diag([-1 1 1]);
    I = repmat(1:numel(hemi), n, 1);  hs = unique(hemi);
    for s = 1:n
        [Q, Rq] = qr(randn(rs, 3));  Q = Q * diag(sign(diag(Rq)));  if det(Q) < 0, Q(:, 1) = -Q(:, 1); end
        for h = hs'
            Qg = Q;  if h == "R", Qg = Fx * Q * Fx; end
            k = find(hemi == h);  I(s, k) = k(knnsearch(U(k, :), U(k, :) * Qg));
        end
    end
end

% Author: Diellor Basha, 2026
