function [T, X, Ex, W] = measure_cptrack(name, opts)
% SCALE.MEASURE_CPTRACK  A drifting rotor and a waxing/waning source, tracked by rheome.detect.track, against the truth.
%
%   [T, X, Ex] = rheome.scale.measure_cptrack(name)           % sphere + the left hemisphere of name
%   [T, X, Ex] = rheome.scale.measure_cptrack('')             % the sphere alone
%   [T, X, Ex] = rheome.scale.measure_cptrack(name, SpeedMS=0.05, Hemis=["L" "R"])
%
% MS1 capability 2.3 (critical points and their curl and divergence over time), the figure panel. Two planted
% critical points on one surface, no instrument (no gain, inverse or noise):
%   rotor   J_r = n x grad g,  g a Gaussian bump (SigmaMM, graph geodesic) whose centre drifts along a
%           geodesic path at SpeedMS: a vortex, curl . n = Lap g at its core, div 0
%   source  J_s = -a(t) grad h, h a fixed bump at the vertex farthest from the path, a(t) = sin^2 rising
%           and falling over the middle half of the record: a source that strengthens then fades,
%           div = -a Lap h at its centre, curl 0
% J = J_r + J_s at Fps for NFrames, linked by rheome.detect.track (types vortex + source; one threshold for
% its strong and weak bands, MinFrac x the rotor's true peak |curl|, so the source's true on/off frames are
% where its true |div| crosses it). The truth for curl and div is closed form: the Laplacian of a geodesic
% Gaussian exp(-r^2/2s^2) at its centre is -2/s^2 on any smooth surface, so curl = -2/s^2 at the rotor core
% and div = +2 a(t)/s^2 at the source centre (both bumps peak at 1). It is independent of the mesh and of
% the face-gradient curl and divergence the detector reads.
%
% X, one row per surface: surface, hemi, speed_ms (planted, geodesic), speed_true_euclid_ms (the core's
%   chord speed, what detect.track's positions measure), rotor_cover (fraction of frames the rotor track
%   holds), rotor_err_mm_median / _max (track position vs the true core), rotor_speednet_ms (track net
%   speed), rotor_curl_relerr (curl at the tracked core vs truth over the tracked frames), rotor_div_leak
%   (|div| / |curl| there; truth 0), source_on_true / source_off_true (frames a(t) x peak >= MinStrength),
%   source_on / source_off (the source track's birth and death frames), source_err_mm_median,
%   source_div_relerr, source_curl_leak, n_other_tracks (confirmed tracks that are neither).
%   Like-for-like speeds over the rotor track's own frames f (lifetime T = f(end) - f(1)):
%   rotor_speedpath_ms (the track's pathLength / T) vs speed_true_path_ms (summed chords of the true core
%   over f / T), and rotor_speednet_ms (net chord / T) vs speed_true_net_ms (|P(f(end)) - P(f(1))| / T).
%   speed_true_euclid_ms (the whole record's summed chords) is kept for continuity only: it is a path speed
%   and must not be read against the net speed (on a folded cortex the 60 mm path bends).
%   rotor_curl_relerr_mesh / source_div_relerr_mesh: the same reads against the cotan Laplacian of the
%   field actually planted (curl n.(n x grad g) = Lap g, div(-a grad h) = -a Lap h), which separates the
%   operator's error from the plant's (edge-graph distance is longer than the geodesic, so the closed form
%   -2/s^2 is not the planted bump's Laplacian on an irregular mesh).
% W: the scale sweep (SweepSigmaMM, static, the rotor at its first true core and the source at q), one row
%   per surface x plant x sigma; plant "edge" (rheome.geom.edgegraph, as the track run) or "heat" (exact
%   geodesic, rheome.geom.geodesic): curl_relerr (vs -2/s^2), curl_relerr_mesh (vs the cotan Lap g),
%   plant_relerr (cotan Lap g vs -2/s^2: the plant's own departure), the same three for the source's div,
%   div_leak (|div| / |curl| at the core), curl_leak (|curl| / |div| at the source), and the whole field
%   against the cotan Laplacian, area-weighted: curl_field_relerr, div_field_relerr, and the same with both
%   sides low-passed to wavelengths > LowpassMM on the cortex's Laplace-Beltrami modes (_lp; NaN on the sphere).
%   curl_weak_relerr / div_weak_relerr: the bump-weighted integral <g, curl> against <g, Lap g> = -int |grad g|^2
%   (the planted field's Dirichlet energy, exact for the P1 field and -pi for a smooth geodesic Gaussian at any
%   sigma): the operator read at the bump's own scale, free of the plant's vertex-scale roughness.
%   rotor_curl_relerr_lp / source_div_relerr_lp in X: the track's reads with both sides low-passed so.
% T: rheome.scale.rows("cptrack"), band = the surface. Ex: per surface the truth and tracked paths, the
%   curl / div series and the source amplitude, for the figure, and .tracks (every confirmed track: type,
%   birth, death, n_frames, median distance to the rotor's true core and to the source, mm).
%
% ⚠ The planted speed is geodesic; detect.track's positions are 3-D (chord), so the truths it is read
%   against are chord-based too (speed_true_*). Path speed carries position jitter upward (detect.track);
%   net speed does not but shrinks with the path's bending: read each against its own truth.
%
% See also: rheome.detect.track, rheome.detect.criticalPoints, rheome.scale.measure_ambient
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        opts.SpeedMS (1,1) double {mustBePositive} = 0.03
        opts.SigmaMM (1,1) double {mustBePositive} = 15
        opts.NFrames (1,1) double {mustBeInteger} = 41
        opts.Fps (1,1) double {mustBePositive} = 20
        opts.MinFrac (1,1) double = 0.5
        opts.Hemis string = "L"
        opts.SphereRadiusMM (1,1) double {mustBePositive} = 100
        opts.Seed (1,1) double = 11
        opts.SweepSigmaMM double = [10 15 20 30 40]
        opts.LowpassMM (1,1) double {mustBePositive} = 60
    end
    rs = rng(opts.Seed);  c = onCleanup(@() rng(rs));
    R = opts.SphereRadiusMM / 1e3;
    [V0, F0] = rheome.geom.icosphere(5);  V0 = V0 ./ vecnorm(V0, 2, 2);
    surf = struct('surface', "sphere", 'hemi', "", 'S', struct('Vertices', R*V0, 'Faces', F0), 'lbo', []);
    if ~isempty(name)
        B = rheome.load.bases(name);
        for h = opts.Hemis(:)', surf(end+1) = struct('surface', "cortex", 'hemi', h, 'S', B.(h).S, 'lbo', B.(h).lbo); end %#ok<AGROW>
    end
    nT = opts.NFrames;  dt = 1 / opts.Fps;  tt = (0:nT-1)' * dt;  sg = opts.SigmaMM / 1e3;
    X = table();  W = table();  T = rheome.scale.rows("cptrack", strings(0,1), [], "");  Ex = struct();
    for s = surf
        V = double(s.S.Vertices);  F = double(s.S.Faces);  nV = size(V, 1);
        Sx = struct('Vertices', V, 'Faces', F, 'nV', nV);
        fg = rheome.operators.face_gradient(V, F);
        n = fg.W * fg.FaceNormal;  n = n ./ vecnorm(n, 2, 2);
        grd = @(x) i_tan([fg.W*(fg.Gx*x) fg.W*(fg.Gy*x) fg.W*(fg.Gz*x)], n);
        G = rheome.geom.edgegraph(Sx);
        [Lc, Mc] = rheome.operators.laplace_beltrami(V, F, 'lumped');  ac = full(diag(Mc));
        lap = @(x) -(Lc * x) ./ ac;                                    % geometric Laplacian (Lc is the PSD stiffness)
        lp = @(x) NaN(size(x));                                        % keep wavelengths > LowpassMM (as measure_ambient)
        if ~isempty(s.lbo)
            k = s.lbo.Lambda(:) <= (2*pi / (opts.LowpassMM/1e3))^2;  Pl = s.lbo.Phi(:, k);
            lp = @(x) Pl * (Pl' * (s.lbo.Mass * x));
        end
        % the rotor's path: a geodesic from a random vertex, as long as the record needs
        cvx = rheome.geom.cortexmask(s.S);  cvi = find(cvx);  len = opts.SpeedMS * tt(end);  v0 = cvi(randi(numel(cvi)));  d0 = distances(G, v0)';
        [~, v1] = min(abs(d0 - len) + 1e9 * ~cvx);  pv = shortestpath(G, v0, v1);   % both ends on cortex
        D = distances(G, pv)';  arc = D(pv(:), 1);  arc = arc(:);   % arclength of each path vertex from v0
        % the source: the vertex farthest from the whole path
        [~, q] = max(min(D, [], 2));  h = exp(-distances(G, q)'.^2 / (2*sg^2));  gh = grd(h);
        tOn = 0.25 * tt(end);  tOff = 0.75 * tt(end);
        amp = (tt >= tOn & tt <= tOff) .* sin(pi * (tt - tOn) / (tOff - tOn)).^2;
        J = zeros(3*nV, nT);  P = zeros(nT, 3);  curlTrue = zeros(nT, 1);  cv = zeros(nT, 1);  LG = zeros(nV, nT);
        for k = 1:nT
            l = min(opts.SpeedMS * tt(k), arc(end));  i = min(find(arc <= l, 1, 'last'), numel(pv) - 1);
            w = (l - arc(i)) / max(arc(i+1) - arc(i), eps);
            dk = (1 - w) * D(:, i) + w * D(:, i+1);  g = exp(-dk.^2 / (2*sg^2));
            P(k, :) = (1 - w) * V(pv(i), :) + w * V(pv(i+1), :);
            [~, cv(k)] = min(vecnorm(V - P(k, :), 2, 2));
            J3 = cross(n, grd(g), 2) - amp(k) * gh;  LG(:, k) = lap(g);
            J(:, k) = reshape(J3.', [], 1);
        end
        curlTrue(:) = -2 / sg^2;  divTrue = amp * 2 / sg^2;          % Lap of exp(-r^2/2s^2) at r = 0 is -2/s^2
        minS = opts.MinFrac * 2 / sg^2;                                % one threshold: strong = weak band
        out = rheome.detect.track(J, Sx, struct('types', {{'vortex', 'source'}}, 'dt', dt, 'minStrength', minS, ...
                                                'highStrength', minS));
        % ---- read back ----
        curlM = rheome.differential.curl(J, Sx, fg);  divM = rheome.differential.divergence(J, Sx, fg);
        tr = out.tracks;  fr = 1:nT;
        [ir, er] = i_pick(tr, "vortex", @(f) P(f, :));
        [is, es] = i_pick(tr, "source", @(f) repmat(V(q, :), numel(f), 1));
        on = fr(abs(divTrue) >= minS);  onT = [min([on NaN]) max([on NaN])];
        r = nan(1, 6);  sx = nan(1, 6);  rtr = struct('frames', [], 'pos', zeros(0, 3));  stv = rtr;  sp = nan(1, 7);  lh = lap(h);
        if ~isnan(ir)
            rtr = tr(ir);  f = rtr.frames(:);  vk = i_near(V, rtr.pos);
            cm = curlM(sub2ind(size(curlM), vk, f));  dm = divM(sub2ind(size(divM), vk, f));
            r = [numel(f) / nT, median(er) * 1e3, max(er) * 1e3, rtr.speedNet, ...
                 norm(cm - curlTrue(f)) / norm(curlTrue(f)), norm(dm) / norm(cm)];
            Tl = (f(end) - f(1)) * dt;  cr = LG(sub2ind(size(LG), vk, f));
            sp = [sum(vecnorm(diff(P(f, :)), 2, 2)) / Tl, rtr.pathLength / Tl, norm(P(f(end), :) - P(f(1), :)) / Tl, ...
                  norm(cm - cr) / norm(cr), NaN, NaN, NaN];
            cl = lp(curlM);  rl = lp(LG);  i2 = sub2ind(size(cl), vk, f);
            sp(6) = norm(cl(i2) - rl(i2)) / norm(rl(i2));
        end
        if ~isnan(is)
            stv = tr(is);  f = stv.frames(:);  dm = divM(q, f)';  cm = curlM(q, f)';
            sx = [stv.birthFrame, stv.deathFrame, median(es) * 1e3, norm(dm - divTrue(f)) / norm(divTrue(f)), ...
                  norm(cm) / norm(dm), NaN];
            sp(5) = norm(dm + amp(f) * lh(q)) / norm(amp(f) * lh(q));
            dl = lp(divM(:, f));  hl = lp(lh);  sp(7) = norm(dl(q, :)' + amp(f) * hl(q)) / norm(amp(f) * hl(q));
        end
        vEu = sum(vecnorm(diff(P), 2, 2)) / tt(end);
        nOther = numel(tr) - ~isnan(ir) - ~isnan(is);
        row = {s.surface, s.hemi, opts.SpeedMS, vEu, r(1), r(2), r(3), r(4), r(5), r(6), ...
               onT(1), onT(2), sx(1), sx(2), sx(3), sx(4), sx(5), nOther, sp(1), sp(2), sp(3), sp(4), sp(5), sp(6), sp(7)};
        X = [X; row]; %#ok<AGROW>
        b = s.surface + iif(s.hemi == "", "", "-" + s.hemi);
        T = [T; rheome.scale.rows("cptrack", ["rotor_cover" "rotor_err_mm_median" "rotor_speednet_ms" ...
             "speed_true_euclid_ms" "rotor_curl_relerr" "rotor_div_leak" "source_on_err_frames" ...
             "source_off_err_frames" "source_err_mm_median" "source_div_relerr" "source_curl_leak" "n_other_tracks" ...
             "speed_true_path_ms" "rotor_speedpath_ms" "speed_true_net_ms" "rotor_curl_relerr_mesh" "source_div_relerr_mesh" ...
             "rotor_curl_relerr_lp" "source_div_relerr_lp"], ...
             [r(1) r(2) r(4) vEu r(5) r(6) sx(1)-onT(1) sx(2)-onT(2) sx(3) sx(4) sx(5) nOther sp], ...
             ["fraction" "mm" "m/s" "m/s" "ratio" "ratio" "frames" "frames" "mm" "ratio" "ratio" "tracks" ...
              "m/s" "m/s" "m/s" "ratio" "ratio" "ratio" "ratio"], b)]; %#ok<AGROW>
        Ws = i_sweep(Sx, fg, n, grd, G, lap, lp, ac, cv(1), q, opts.SweepSigmaMM / 1e3);
        Ws = [table(repmat(s.surface, height(Ws), 1), repmat(s.hemi, height(Ws), 1), 'VariableNames', {'surface', 'hemi'}) Ws];
        W = [W; Ws]; %#ok<AGROW>
        for i = 1:height(Ws)
            T = [T; rheome.scale.rows("cptrack", ["sweep_curl_relerr" "sweep_curl_relerr_mesh" "sweep_curl_field_relerr_lp" ...
                 "sweep_div_relerr" "sweep_div_relerr_mesh" "sweep_div_field_relerr_lp"], [Ws.curl_relerr(i) Ws.curl_relerr_mesh(i) ...
                 Ws.curl_field_relerr_lp(i) Ws.div_relerr(i) Ws.div_relerr_mesh(i) Ws.div_field_relerr_lp(i)], "ratio", ...
                 b + "/" + Ws.plant(i) + "/" + Ws.sigma_mm(i) + " mm")]; %#ok<AGROW>
        end
        Ttr = i_tracks(tr, P, V(q, :));
        Ex.(matlab.lang.makeValidName(char(s.surface + s.hemi))) = struct('surface', char(s.surface), 'hemi', char(s.hemi), ...
            'vertices_mm', 1e3*V, 'faces', F, 't_s', tt, 'core_true_mm', 1e3*P, 'core_vertex_true', cv, ...
            'rotor_frames', rtr.frames(:), 'rotor_mm', 1e3*rtr.pos, 'source_vertex', q, 'source_amp', amp, ...
            'source_frames', stv.frames(:), 'source_mm', 1e3*stv.pos, 'curl_true', curlTrue, 'div_true', divTrue, ...
            'curl_at_true_core', curlM(sub2ind(size(curlM), cv, fr')), 'div_at_source', divM(q, :)', ...
            'min_strength', minS, 'tracks', Ttr, 'units', 'positions mm; curl, div 1/m x J (J arbitrary)');
    end
    X.Properties.VariableNames = {'surface','hemi','speed_ms','speed_true_euclid_ms','rotor_cover', ...
        'rotor_err_mm_median','rotor_err_mm_max','rotor_speednet_ms','rotor_curl_relerr','rotor_div_leak', ...
        'source_on_true','source_off_true','source_on','source_off','source_err_mm_median','source_div_relerr', ...
        'source_curl_leak','n_other_tracks','speed_true_path_ms','rotor_speedpath_ms','speed_true_net_ms', ...
        'rotor_curl_relerr_mesh','source_div_relerr_mesh','rotor_curl_relerr_lp','source_div_relerr_lp'};
end

function W = i_sweep(Sx, fg, n, grd, G, lap, lp, a, c, q, sigmas)
% the operators at a scale: a static rotor at c and source at q, planted with edge-graph or exact geodesic distance
    emb = @(w) reshape(w.', [], 1);  W = table();  fe = @(x, y) sqrt(a' * (x - y).^2) / sqrt(a' * y.^2);   % area-weighted
    [dh, sol] = rheome.geom.geodesic(Sx, c);  dq = rheome.geom.geodesic(Sx, q, sol);
    D = struct('edge', {{distances(G, c)', distances(G, q)'}}, 'heat', {{dh, dq}});
    for plant = ["edge" "heat"]
        d = D.(plant);
        for sg = sigmas(:)'
            g = exp(-d{1}.^2 / (2*sg^2));  h = exp(-d{2}.^2 / (2*sg^2));
            Jr = emb(cross(n, grd(g), 2));  Js = emb(-grd(h));
            cr = rheome.differential.curl(Jr, Sx, fg);  dr = rheome.differential.divergence(Jr, Sx, fg);
            cs = rheome.differential.curl(Js, Sx, fg);  ds = rheome.differential.divergence(Js, Sx, fg);
            lg = lap(g);  lh = lap(h);  ref = 2 / sg^2;
            W = [W; table(plant, 1e3*sg, abs(cr(c) + ref) / ref, abs(cr(c) - lg(c)) / abs(lg(c)), abs(lg(c) + ref) / ref, ...
                 abs(ds(q) - ref) / ref, abs(ds(q) + lh(q)) / abs(lh(q)), abs(-lh(q) - ref) / ref, ...
                 abs(dr(c)) / abs(cr(c)), abs(cs(q)) / abs(ds(q)), fe(cr, lg), fe(lp(cr), lp(lg)), fe(ds, -lh), fe(lp(ds), lp(-lh)), ...
                 abs(a' * (g .* (cr - lg))) / abs(a' * (g .* lg)), abs(a' * (h .* (ds + lh))) / abs(a' * (h .* lh)), ...
                 'VariableNames', {'plant', 'sigma_mm', 'curl_relerr', 'curl_relerr_mesh', 'curl_plant_relerr', 'div_relerr', ...
                 'div_relerr_mesh', 'div_plant_relerr', 'div_leak', 'curl_leak', 'curl_field_relerr', 'curl_field_relerr_lp', ...
                 'div_field_relerr', 'div_field_relerr_lp', 'curl_weak_relerr', 'div_weak_relerr'})]; %#ok<AGROW>
        end
    end
end

function Tt = i_tracks(tr, P, pq)
% every confirmed track and how far it sits from each planted critical point (median over its frames, mm)
    Tt = table('Size', [numel(tr) 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'type', 'birth', 'death', 'n_frames', 'rotor_dist_mm', 'source_dist_mm'});
    for i = 1:numel(tr)
        f = tr(i).frames(:);
        Tt(i, :) = {string(tr(i).type), f(1), f(end), numel(f), 1e3*median(vecnorm(tr(i).pos - P(f, :), 2, 2)), ...
                    1e3*median(vecnorm(tr(i).pos - pq, 2, 2))};
    end
end

function [best, err] = i_pick(tr, type, truth)
% the track of this type nearest its truth (median distance over its frames), or NaN
    best = NaN;  err = [];  m = Inf;
    for i = 1:numel(tr)
        if ~strcmp(tr(i).type, type), continue; end
        e = vecnorm(tr(i).pos - truth(tr(i).frames(:)), 2, 2);
        if median(e) < m, m = median(e);  best = i;  err = e; end
    end
end

function k = i_near(V, p)
    k = zeros(size(p, 1), 1);
    for i = 1:size(p, 1), [~, k(i)] = min(vecnorm(V - p(i, :), 2, 2)); end
end

function g = i_tan(g, n)
    g = g - sum(g .* n, 2) .* n;
end

function z = iif(c, a, b)
    if c, z = a; else, z = b; end
end

% Author: Diellor Basha, 2026
