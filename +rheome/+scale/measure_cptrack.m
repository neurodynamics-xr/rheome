function [T, X, Ex] = measure_cptrack(name, opts)
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
% T: rheome.scale.rows("cptrack"), band = the surface. Ex: per surface the truth and tracked paths, the
%   curl / div series and the source amplitude, for the figure.
%
% ⚠ The path speed is geodesic; detect.track's speed is from 3-D positions (chord). They differ by the
%   folding, so both are reported (see rheome.flow.movingvortex).
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
    end
    rs = rng(opts.Seed);  c = onCleanup(@() rng(rs));
    R = opts.SphereRadiusMM / 1e3;
    [V0, F0] = rheome.geom.icosphere(5);  V0 = V0 ./ vecnorm(V0, 2, 2);
    surf = struct('surface', "sphere", 'hemi', "", 'S', struct('Vertices', R*V0, 'Faces', F0));
    if ~isempty(name)
        B = rheome.load.bases(name);
        for h = opts.Hemis(:)', surf(end+1) = struct('surface', "cortex", 'hemi', h, 'S', B.(h).S); end %#ok<AGROW>
    end
    nT = opts.NFrames;  dt = 1 / opts.Fps;  tt = (0:nT-1)' * dt;  sg = opts.SigmaMM / 1e3;
    X = table();  T = rheome.scale.rows("cptrack", strings(0,1), [], "");  Ex = struct();
    for s = surf
        V = double(s.S.Vertices);  F = double(s.S.Faces);  nV = size(V, 1);
        Sx = struct('Vertices', V, 'Faces', F, 'nV', nV);
        fg = rheome.operators.face_gradient(V, F);
        n = fg.W * fg.FaceNormal;  n = n ./ vecnorm(n, 2, 2);
        grd = @(x) i_tan([fg.W*(fg.Gx*x) fg.W*(fg.Gy*x) fg.W*(fg.Gz*x)], n);
        G = rheome.geom.edgegraph(Sx);
        % the rotor's path: a geodesic from a random vertex, as long as the record needs
        cvx = rheome.geom.cortexmask(s.S);  cvi = find(cvx);  len = opts.SpeedMS * tt(end);  v0 = cvi(randi(numel(cvi)));  d0 = distances(G, v0)';
        [~, v1] = min(abs(d0 - len) + 1e9 * ~cvx);                    % both ends on cortex  pv = shortestpath(G, v0, v1);
        D = distances(G, pv)';  arc = D(pv(:), 1);  arc = arc(:);   % arclength of each path vertex from v0
        % the source: the vertex farthest from the whole path
        [~, q] = max(min(D, [], 2));  h = exp(-distances(G, q)'.^2 / (2*sg^2));  gh = grd(h);
        tOn = 0.25 * tt(end);  tOff = 0.75 * tt(end);
        amp = (tt >= tOn & tt <= tOff) .* sin(pi * (tt - tOn) / (tOff - tOn)).^2;
        J = zeros(3*nV, nT);  P = zeros(nT, 3);  curlTrue = zeros(nT, 1);  cv = zeros(nT, 1);
        for k = 1:nT
            l = min(opts.SpeedMS * tt(k), arc(end));  i = min(find(arc <= l, 1, 'last'), numel(pv) - 1);
            w = (l - arc(i)) / max(arc(i+1) - arc(i), eps);
            dk = (1 - w) * D(:, i) + w * D(:, i+1);  g = exp(-dk.^2 / (2*sg^2));
            P(k, :) = (1 - w) * V(pv(i), :) + w * V(pv(i+1), :);
            [~, cv(k)] = min(vecnorm(V - P(k, :), 2, 2));
            J3 = cross(n, grd(g), 2) - amp(k) * gh;
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
        r = nan(1, 6);  sx = nan(1, 6);  rtr = struct('frames', [], 'pos', zeros(0, 3));  stv = rtr;
        if ~isnan(ir)
            rtr = tr(ir);  f = rtr.frames(:);  vk = i_near(V, rtr.pos);
            cm = curlM(sub2ind(size(curlM), vk, f));  dm = divM(sub2ind(size(divM), vk, f));
            r = [numel(f) / nT, median(er) * 1e3, max(er) * 1e3, rtr.speedNet, ...
                 norm(cm - curlTrue(f)) / norm(curlTrue(f)), norm(dm) / norm(cm)];
        end
        if ~isnan(is)
            stv = tr(is);  f = stv.frames(:);  dm = divM(q, f)';  cm = curlM(q, f)';
            sx = [stv.birthFrame, stv.deathFrame, median(es) * 1e3, norm(dm - divTrue(f)) / norm(divTrue(f)), ...
                  norm(cm) / norm(dm), NaN];
        end
        vEu = sum(vecnorm(diff(P), 2, 2)) / tt(end);
        nOther = numel(tr) - ~isnan(ir) - ~isnan(is);
        row = {s.surface, s.hemi, opts.SpeedMS, vEu, r(1), r(2), r(3), r(4), r(5), r(6), ...
               onT(1), onT(2), sx(1), sx(2), sx(3), sx(4), sx(5), nOther};
        X = [X; row]; %#ok<AGROW>
        b = s.surface + iif(s.hemi == "", "", "-" + s.hemi);
        T = [T; rheome.scale.rows("cptrack", ["rotor_cover" "rotor_err_mm_median" "rotor_speednet_ms" ...
             "speed_true_euclid_ms" "rotor_curl_relerr" "rotor_div_leak" "source_on_err_frames" ...
             "source_off_err_frames" "source_err_mm_median" "source_div_relerr" "source_curl_leak" "n_other_tracks"], ...
             [r(1) r(2) r(4) vEu r(5) r(6) sx(1)-onT(1) sx(2)-onT(2) sx(3) sx(4) sx(5) nOther], ...
             ["fraction" "mm" "m/s" "m/s" "ratio" "ratio" "frames" "frames" "mm" "ratio" "ratio" "tracks"], b)]; %#ok<AGROW>
        Ex.(matlab.lang.makeValidName(char(s.surface + s.hemi))) = struct('surface', char(s.surface), 'hemi', char(s.hemi), ...
            'vertices_mm', 1e3*V, 'faces', F, 't_s', tt, 'core_true_mm', 1e3*P, 'core_vertex_true', cv, ...
            'rotor_frames', rtr.frames(:), 'rotor_mm', 1e3*rtr.pos, 'source_vertex', q, 'source_amp', amp, ...
            'source_frames', stv.frames(:), 'source_mm', 1e3*stv.pos, 'curl_true', curlTrue, 'div_true', divTrue, ...
            'curl_at_true_core', curlM(sub2ind(size(curlM), cv, fr')), 'div_at_source', divM(q, :)', ...
            'min_strength', minS, 'units', 'positions mm; curl, div 1/m x J (J arbitrary)');
    end
    X.Properties.VariableNames = {'surface','hemi','speed_ms','speed_true_euclid_ms','rotor_cover', ...
        'rotor_err_mm_median','rotor_err_mm_max','rotor_speednet_ms','rotor_curl_relerr','rotor_div_leak', ...
        'source_on_true','source_off_true','source_on','source_off','source_err_mm_median','source_div_relerr', ...
        'source_curl_leak','n_other_tracks'};
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
