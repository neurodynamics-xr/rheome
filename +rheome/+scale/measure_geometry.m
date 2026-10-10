function [T, X, Z] = measure_geometry(name, S, opts)
% SCALE.MEASURE_GEOMETRY  The atlas geometry checks on one subject's cortex: atom-tile overlap, gauge, roll-up,
% and the size of parcels and tiles.
%
%   [T, X] = rheome.scale.measure_geometry(name)
%   [T, X, Z] = rheome.scale.measure_geometry(name)     % Z.parcels, Z.tiles: G12g, G12h
%   [T, X] = rheome.scale.measure_geometry(name, rheome.scale.sensors(name), Levels=2:5, RollupDepths=3:6)
%
% The geometry parts of MS1 group plan G12, per hemisphere, on rheome.geom.tree (area bisection):
%
% 1. WAVELET-TILE MATCH (Fig. 2A). At each level l, one mexhat atom g(x) = x e^-x per tile at its most
%    interior vertex (largest geodesic distance to the tile's boundary), t_l = A_l / (4 pi) with
%    A_l = hemisphere area / 2^l, so the flat-space central lobe (area 4 pi t) equals one tile.
%    psi = Phi (g(t lambda) .* Phi(c,:)'). Per tile: lobe_ratio = area of the connected positive lobe
%    containing c / tile area; energy_in_tile = sum over the tile of a psi^2 / sum of a psi^2 (lumped area a).
%    Rows (band "L<l>"): lobe_ratio and energy_in_tile, medians over the tiles of both hemispheres.
% 2. GAUGE (Fig. 16). rheome.operators.gauge Method "trivial", +1 at the two faces of largest separation
%    (the default). Rows (band "L" | "R"): n_singular (expected 2), charge_sum (expected chi = 2),
%    neighbour_angle_deg (median between neighbouring frame vectors), residual_rad (own connection,
%    0 by construction) and residual_lc_rad (against Levi-Civita).
% 3. ROLL-UP EXACTNESS (Section 13.4). An 8-16 Hz co-power connectome C = Re(sum_xyz X X^H) of the
%    tile-summed minimum-norm current (S.Res on the whole cached recording, analytic band-pass at
%    the sensors), direct at each depth in RollupDepths against block sums of the deepest. Rows (band
%    "L" | "R"): rollup_max_rel_error (max over depths of max|rolled - direct| / max|direct|).
% 4. PARCEL AND TILE SIZE (G12g, G12h). One ruler for both: the equivalent-disc diameter d_eq = 2 sqrt(A/pi),
%    A the lumped vertex area (1/3 of the incident triangles) summed over the members, on this cortex.
%    Z.parcels: one row per parcel of each of opts.Atlases (Brainstorm's scouts; any unknown /
%    corpuscallosum / Medial_wall label dropped) -- atlas, label, hemi, n_vertices, area_mm2, d_eq_mm.
%    Z.tiles: one row per rheome.geom.tree tile at depths opts.SizeDepths (depth 0 = the hemisphere) --
%    hemi, depth, node_id, n_vertices, area_mm2, d_eq_mm, geodesic_diameter_mm (double-sweep on the
%    mesh edge graph: a lower bound on the longest geodesic in the tile, and specific to this Fiedler tiling).
%    Rows: hemisphere_area (cm2, band "L" | "R"); parcel_deq_median, n_parcels (band = atlas, over both
%    hemispheres); tile_deq_median, tile_area_median (band "D<depth>"). A missing atlas gives n_parcels 0.
%    The fraction of parcels below 2 r50 needs r50 (cf-scale) and is taken at the group reduction.
%
% ⚠ THE TILING IS rheome.geom.tree (Fiedler bisection), NOT THE ATLAS's EQUAL-AREA GEODESIC BISECTION that
% drew Fig. 2 (nsp atlas ladder, F13 card: lobe 0.83/0.95/0.97/0.97, energy 80/84/85/86 % at l = 2-5,
% 700 modes). Both halve area per level; tile shapes differ, so the numbers are comparable, not identical.
% The basis is the cache's own (rheome.load.bases: 1000 modes per hemisphere in rheome.scale.importsubject).
% ⚠ The gauge's poles here are geometric (largest separation); MS1 Fig. 16B put them at the two
% lowest-alpha faces. The singularity count and the smoothness do not depend on
% where they are; the canonical (registered-sphere) frame is the atlas's (group plan P8).
%
% See also: rheome.geom.tree, rheome.geom.tiles, rheome.operators.gauge, rheome.ingest.rollup
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Levels double = 2:5
        opts.RollupDepths double = 3:6
        opts.Band (1,2) double = [8 16]
        opts.Hemis string = ["L" "R"]
        opts.Atlases string = ["Desikan-Killiany" "Destrieux"]
        opts.SizeDepths double = 0:3
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', opts.Band(1), ...
                    'HalfPowerFrequency2', opts.Band(2), 'SampleRate', fs);
    Fa = hilbert(filtfilt(bp, F.')).';  clear F
    T = rheome.scale.rows("geometry", strings(0,1), [], "");  X = table();  Zt = table();  ag = [];
    dmax = max([opts.Levels opts.RollupDepths opts.SizeDepths]);
    for h = opts.Hemis(:)'
        Hh = S.B.(h);  Sh = Hh.S;  lbo = Hh.lbo;  V = double(Sh.Vertices);  Fc = double(Sh.Faces);
        a = full(sum(lbo.Mass, 2));  Atot = sum(a);  lam = lbo.Lambda(:);  P = lbo.Phi;
        Tr = rheome.geom.tree(Sh, L=lbo.L, M=lbo.M, MaxDepth=dmax);
        E = unique(sort([Fc(:,[1 2]); Fc(:,[2 3]); Fc(:,[3 1])], 2), 'rows');
        Gm = graph(E(:,1), E(:,2), vecnorm(V(E(:,1),:) - V(E(:,2),:), 2, 2), size(V, 1));
        ag(Hh.gv) = a; %#ok<AGROW>                            % the whole cortex's vertex areas, for the parcels

        % 4. tile size (G12h): the tree's own area and equivalent-disc diameter
        for r = find(ismember(Tr.depth, opts.SizeDepths))'
            in = Tr.members{r}(:);  Gs = subgraph(Gm, in);
            d1 = distances(Gs, 1);  d1(~isfinite(d1)) = -1;  [~, f] = max(d1);   % double sweep from vertex 1
            dg = distances(Gs, f);  dg = max(dg(isfinite(dg)));
            Zt = [Zt; {h, Tr.depth(r), Tr.node_id(r), numel(in), 1e6 * Tr.area(r), 1e3 * Tr.diameter(r), 1e3 * dg}]; %#ok<AGROW>
        end
        T = [T; rheome.scale.rows("geometry", "hemisphere_area", 1e4 * Atot, "cm2", h)]; %#ok<AGROW>

        % 1. wavelet-tile match
        for l = opts.Levels
            G = rheome.geom.tiles(Tr, Sh, l);  t = Atot / 2^l / (4*pi);  g = (lam*t) .* exp(-lam*t);
            bnd = false(size(V,1), 1);  ab = G.tileOf(E(:,1)) ~= G.tileOf(E(:,2));  bnd(E(ab, :)) = true;
            for s = 1:numel(G.node_id)
                in = find(G.tileOf == s);  if numel(in) < 3, continue; end
                Gs = subgraph(Gm, in);  src = find(bnd(in));
                if isempty(src), c = in(1); else, d = min(distances(Gs, src), [], 1); [~, k] = max(d); c = in(k); end
                psi = P * (g .* P(c, :)');
                pos = find(psi > 0);  lobe = [];
                if psi(c) > 0, cc = conncomp(subgraph(Gm, pos));  lobe = pos(cc == cc(pos == c)); end
                X = [X; {h, l, s, c, numel(in), sum(a(lobe)) / sum(a(in)), sum(a(in) .* psi(in).^2) / sum(a .* psi.^2)}]; %#ok<AGROW>
            end
        end

        % 2. gauge
        gT = rheome.operators.gauge(V, Fc, Method="trivial");
        T = [T; rheome.scale.rows("geometry", ["n_singular" "charge_sum" "neighbour_angle_deg" "residual_rad" "residual_lc_rad"], ...
             [numel(gT.singular) sum(gT.charge) gT.neighbourAngle gT.residual gT.residualLC], ...
             ["faces" "index" "deg" "rad" "rad"], h)]; %#ok<AGROW>

        % 3. roll-up exactness of the tile co-power connectome
        gv = double(Hh.gv(:));  dd = max(opts.RollupDepths);  Gd = rheome.geom.tiles(Tr, Sh, dd);
        if ~Gd.heap, error('scale:geometry:heap', 'hemisphere %s: tree not full to depth %d, ancestors are not arithmetic', h, dd); end
        Xd = 0;  Cd = 0;
        for k = 1:3                                            % x, y, z: tile sums of each current component
            Xk = (double(Gd.P)' * S.Res.ImagingKernel((gv - 1)*3 + k, :)) * Fa;
            Cd = Cd + real(Xk * Xk');
        end
        err = 0;
        for l = setdiff(opts.RollupDepths, dd)
            Gl = rheome.geom.tiles(Tr, Sh, l);  [~, par] = ismember(floor(Gd.node_id / 2^(dd - l)), Gl.node_id);
            Ag = sparse(1:numel(par), par, 1, numel(par), numel(Gl.node_id));
            Cl = 0;
            for k = 1:3
                Xk = (double(Gl.P)' * S.Res.ImagingKernel((gv - 1)*3 + k, :)) * Fa;  Cl = Cl + real(Xk * Xk');
            end
            err = max(err, max(abs(Ag' * Cd * Ag - Cl), [], 'all') / max(abs(Cl), [], 'all'));
        end
        T = [T; rheome.scale.rows("geometry", "rollup_max_rel_error", err, "relative", h)]; %#ok<AGROW>
        fprintf('[geometry %s] %s: %d singular faces, %.1f deg, roll-up %.2e\n', name, h, numel(gT.singular), gT.neighbourAngle, err);
    end
    X.Properties.VariableNames = {'hemi','level','tile','centre','n_vertices','lobe_ratio','energy_in_tile'};
    Zt.Properties.VariableNames = {'hemi','depth','node_id','n_vertices','area_mm2','d_eq_mm','geodesic_diameter_mm'};
    for d = opts.SizeDepths
        k = Zt.depth == d;
        T = [T; rheome.scale.rows("geometry", ["tile_deq_median" "tile_area_median"], ...
             [median(Zt.d_eq_mm(k)) median(Zt.area_mm2(k))], ["mm" "mm2"], "D" + d)]; %#ok<AGROW>
    end

    % 4. parcel size (G12g)
    Zp = table('Size', [0 6], 'VariableTypes', {'string','string','string','double','double','double'}, ...
               'VariableNames', {'atlas','label','hemi','n_vertices','area_mm2','d_eq_mm'});
    for at = opts.Atlases(:)'
        try, A = rheome.load.atlas(name, at);
        catch e, warning('scale:geometry:atlas', '%s: %s', at, e.message);  A = struct('Label', {{}}, 'Membership', sparse(0, numel(ag)));
        end
        lab = string(A.Label(:));  keep = ~contains(lower(lab), ["unknown" "corpuscallosum" "medial_wall"]);
        Pm = double(A.Membership(keep, :));  lab = lab(keep);  av = zeros(size(Pm, 2), 1);  av(1:numel(ag)) = ag;
        ar = 1e6 * (Pm * av);  hm = extractAfter(lab, strlength(lab) - 1);  hm(~ismember(hm, ["L" "R"])) = "";
        Zp = [Zp; table(repmat(at, numel(lab), 1), lab, hm, full(sum(Pm, 2)), ar, 2 * sqrt(ar / pi), ...
                        'VariableNames', Zp.Properties.VariableNames)]; %#ok<AGROW>
        T = [T; rheome.scale.rows("geometry", ["parcel_deq_median" "n_parcels"], ...
             [median(2 * sqrt(ar / pi)) numel(lab)], ["mm" "parcels"], at)]; %#ok<AGROW>
    end
    Z = struct('parcels', Zp, 'tiles', Zt);
    for l = opts.Levels
        k = X.level == l;
        T = [T; rheome.scale.rows("geometry", ["lobe_ratio" "energy_in_tile" "n_tiles"], ...
             [median(X.lobe_ratio(k)) median(X.energy_in_tile(k)) nnz(k)], ["ratio" "fraction" "tiles"], ...
             "L" + l)]; %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
