function Op = operator(S, Tp, opts)
% MODES.OPERATOR  Sensors -> group tiles: the normal current, low-passed at a declared wavelength, per tile.
%
%   Op = rheome.modes.operator(S, rheome.scale.template(4))                  % >= 100 mm, depth 3
%   Op = rheome.modes.operator(S, Tp, CutoffMM=130, Depth=2)
%
% Three linear maps in one matrix, so a tile series is W * recording:
%   1. the unconstrained kernel (S.Res, the paper's whitened minimum norm) projected on each vertex's normal:
%      the signed cortical current q, which keeps its sign across a sulcus where the 3-vector flips;
%   2. a scalar Laplace-Beltrami low-pass keeping wavelengths >= CutoffMM: the gain of the TIGHT logitersine
%      bank's members at or above the cutoff, h = sqrt(sum_kept g_m^2) -- the gain rheome.flow.bandlimit uses
%      (h <= 1, h = 1 on the kept band), here on q instead of on the tangential potentials. ⚠ flow.bandlimit
%      itself rebuilds only the TANGENTIAL current (grad Phi + n x grad Psi) and drops the normal part, i.e.
%      the cortical dipole; that is right for flow, wrong for a mode of the current;
%   3. the area-weighted mean over each group tile of depth Depth (rheome.scale.grouptiles: the template's tiles
%      carried through the registered sphere, so tile k is the same place in every participant).
%
% Op.W [nTiles x nSel]   Op.tiles table: hemi, tile (code), area_mm2, x_mm y_mm z_mm (participant centroid),
%   u1 u2 u3 (template sphere centre), tx ty tz (template coordinates, mm; SCS: x anterior, y left, z up)
% Op.A [nTiles x nTiles] adjacency (within hemisphere)
% Op.R [nTiles x nTiles] tile resolution: column i = what W reads from a uniform normal current over tile i,
%   through the participant's gain (S.G) -- the instrument's view of one tile generator (leakage, N3)
% Op.C1 [nSel x nSel] sensor covariance of independent normal sources at every vertex, variance ~ area (N1)
% Op.keptWavelengthMM, Op.cutoffMM, Op.depth
%
% See also: rheome.modes.fields, rheome.flow.bandlimit, rheome.scale.grouptiles
%
% Author: Diellor Basha, 2026

    arguments
        S (1,1) struct
        Tp (1,1) struct
        opts.CutoffMM (1,1) double {mustBePositive} = 100
        opts.Depth (1,1) double {mustBeInteger, mustBePositive} = 3
        opts.Voices (1,1) double {mustBePositive} = 2
        opts.Hemis string = ["L" "R"]
    end
    K = S.Res.ImagingKernel;  W = [];  R = [];  A = [];  tiles = table();  Pall = [];  kept = [];
    Gn = zeros(size(S.G, 1), 0);  C1 = 0;
    for h = opts.Hemis(:)'
        Hh = S.B.(h);  gv = double(Hh.gv(:));  lbo = Hh.lbo;
        nrm = i_normals(S.Sf, gv, Hh.S);
        Kq = nrm(:,1) .* K((gv-1)*3 + 1, :) + nrm(:,2) .* K((gv-1)*3 + 2, :) + nrm(:,3) .* K((gv-1)*3 + 3, :);
        lam = lbo.Lambda(:);
        gfb = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', opts.Voices);
        G = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:gfb.NumMembers, 'uni', 0));
        wl = 1e3 * wavelengths(gfb);  keep = wl >= opts.CutoffMM;
        if ~any(keep), error('modes:operator:cutoff', 'no bank member has a wavelength >= %g mm (coarsest %.0f mm)', opts.CutoffMM, max(wl)); end
        hg = sqrt(sum(G(:, keep).^2, 2));  kept = union(kept, wl(keep));
        [P, Ah, ids, C] = rheome.scale.grouptiles(S, Tp, h, opts.Depth);
        a = full(sum(lbo.Mass, 2));  Pa = (P .* a) ./ full(sum(P .* a, 1));      % area-weighted mean
        Lw = lbo.Phi * (hg .* (lbo.Phi' * (lbo.Mass * Kq)));                    % low-passed q per vertex
        W = [W; full(Pa' * Lw)]; %#ok<AGROW>
        Gh = nrm(:,1)' .* S.G(:, (gv-1)*3 + 1) + nrm(:,2)' .* S.G(:, (gv-1)*3 + 2) + nrm(:,3)' .* S.G(:, (gv-1)*3 + 3);
        Gn = [Gn, Gh * (P .* a)]; %#ok<AGROW>                                    % uniform current density per tile
        C1 = C1 + (Gh .* a') * Gh';                                              % independent sources, variance ~ area
        A = blkdiag(A, full(Ah));
        V = double(Hh.S.Vertices);  xc = full(Pa' * V) * 1e3;
        tiles = [tiles; table(repmat(h, numel(ids), 1), ids(:), full(sum(P .* a, 1))' * 1e6, xc(:,1), xc(:,2), xc(:,3), ...
                 C.u(:,1), C.u(:,2), C.u(:,3), 1e3*C.x(:,1), 1e3*C.x(:,2), 1e3*C.x(:,3), 'VariableNames', ...
                 {'hemi','tile','area_mm2','x_mm','y_mm','z_mm','u1','u2','u3','tx','ty','tz'})]; %#ok<AGROW>
    end
    R = W * Gn;  R = R ./ max(abs(R), [], 1);
    Op = struct('W', W, 'A', logical(A), 'R', R, 'C1', (C1 + C1') / 2, 'tiles', tiles, 'keptWavelengthMM', kept(:)', ...
                'cutoffMM', opts.CutoffMM, 'depth', opts.Depth);
end

function N = i_normals(Sf, gv, Sh)
    if isfield(Sf, 'VertNormals') && ~isempty(Sf.VertNormals), N = double(Sf.VertNormals(gv, :));
    else
        V = double(Sh.Vertices);  Fc = double(Sh.Faces);  N = zeros(size(V));
        fn = cross(V(Fc(:,2),:) - V(Fc(:,1),:), V(Fc(:,3),:) - V(Fc(:,1),:), 2);
        for j = 1:3, N = N + accumarray([repmat(Fc(:,j), 3, 1) kron((1:3)', ones(size(Fc,1), 1))], fn(:), size(V)); end
    end
    N = N ./ vecnorm(N, 2, 2);
end

% Author: Diellor Basha, 2026
