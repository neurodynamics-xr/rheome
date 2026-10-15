function [T, X, Ex] = measure_ambient(name, opts)
% SCALE.MEASURE_AMBIENT  A planted tilted field read ambient vs tangential, on the sphere and the subject's cortex.
%
%   [T, X, Ex] = rheome.scale.measure_ambient(name)            % sphere + the left hemisphere of name
%   [T, X, Ex] = rheome.scale.measure_ambient('')              % the sphere alone (closed form only)
%   [T, X, Ex] = rheome.scale.measure_ambient(name, TiltsDeg=[0 45 90], Hemis=["L" "R"])
%
% MS1 capability 1.4 (ambient-vector fields), the figure panel. The current is a moment, not a tangent
% vector: a patch whose moment tilts away from the normal has a normal share, and on a curved sheet that
% share reads as divergence. The two readings differ by the mean-curvature coupling,
%
%     div J  =  div_S J_t  +  2H (J . n)             (n outward, H = 1/R on a sphere of radius R)
%
% (rheome.operators.weak_differential, about_sensor_to_cortex). This plants J = f (cos(theta) n + sin(theta) t),
% f a Gaussian bump (SigmaMM, graph geodesic, rheome.geom.edgegraph) and t the tangential projection of a fixed
% axis, at each tilt theta from the normal, and reads:
%   D_amb = rheome.differential.divergence(J)       the ambient reading
%   D_tan = rheome.differential.divergence(J_t)     the tangential reading, J_t = J - (J.n) n
%   C     = D_amb - D_tan                           what the normal share adds (exactly div of (J.n) n)
%   K     = 2H (J.n)                                its closed form: H = 1/R on the sphere; on the cortex
%                                                   H from the cotan mean-curvature normal, (M\L x).n / 2
% No instrument: no gain, no inverse, no noise. Geometry and operators only.
%
% X, one row per surface x tilt: surface, hemi, tilt_deg, sigma_mm, centre_vertex, n_vertices,
%   H_relerr (sphere: cotan H vs 1/R), couple_relerr (||C - K|| / ||K||, area-weighted), couple_r
%   (Pearson r of C and K), couple_relerr_60mm / couple_r_60mm (both low-passed to wavelengths > 60 mm on
%   the cortex's Laplace-Beltrami modes; NaN on the sphere), normal_share (||C|| / ||D_amb||: the fraction
%   of the ambient divergence the normal share makes), curl_normal_rel (||curl of the normal part|| / ||C||;
%   closed form 0).
% T: rheome.scale.rows("ambient"), band "<surface>/<tilt> deg": couple_relerr, couple_r, normal_share and
%   (cortex) couple_relerr_60mm. Ex: figure data per surface at ExampleTiltDeg (vertices_mm, faces,
%   normals, H, f, J [nV x 3], D_amb, D_tan, C, K).
%
% ⚠ On the cortex the closed form uses a DIFFERENT discretisation of H (cotan) from the one the divergence
%   carries (hat gradients x vertex normals), so the per-vertex error there is the mesh's, not the
%   identity's: read it at a scale (the 60 mm columns). The sphere is the exact check (tests/tScaleAmbient).
%
% See also: rheome.differential.divergence, rheome.operators.weak_differential, rheome.scale.measure_cptrack
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        opts.TiltsDeg double = [0 30 45 60 90]
        opts.SigmaMM (1,1) double {mustBePositive} = 20
        opts.Hemis string = "L"
        opts.SphereRadiusMM (1,1) double {mustBePositive} = 100
        opts.ExampleTiltDeg (1,1) double = 45
        opts.LowpassMM (1,1) double {mustBePositive} = 60
        opts.Seed (1,1) double = 11
    end
    rs = rng(opts.Seed);  c = onCleanup(@() rng(rs));
    R = opts.SphereRadiusMM / 1e3;
    [V0, F0] = rheome.geom.icosphere(5);  V0 = V0 ./ vecnorm(V0, 2, 2);
    surf = struct('surface', "sphere", 'hemi', "", 'S', struct('Vertices', R*V0, 'Faces', F0), 'lbo', []);
    if ~isempty(name)
        B = rheome.load.bases(name);
        for h = opts.Hemis(:)'
            surf(end+1) = struct('surface', "cortex", 'hemi', h, 'S', B.(h).S, 'lbo', B.(h).lbo); %#ok<AGROW>
        end
    end
    X = table();  T = rheome.scale.rows("ambient", strings(0,1), [], "");  Ex = struct();
    for s = surf
        V = double(s.S.Vertices);  F = double(s.S.Faces);  nV = size(V, 1);  Sx = struct('Vertices', V, 'Faces', F);
        fg = rheome.operators.face_gradient(V, F);
        n = fg.W * fg.FaceNormal;  n = n ./ vecnorm(n, 2, 2);          % the normals the face divergence sees
        [L, M] = rheome.operators.laplace_beltrami(V, F, 'lumped');  a = full(diag(M));
        H = 0.5 * sum(((L * V) ./ a) .* n, 2);                         % cotan mean curvature, 1/m
        if s.surface == "sphere"                                       % H = 1/R, signed as n is (the winding)
            Hk = sign(median(H)) * ones(nV, 1) / R;  Hrel = sqrt(a' * (H - Hk).^2 / sum(a)) * R;
        else, Hk = H;  Hrel = NaN; end
        cvx = find(rheome.geom.cortexmask(s.S));  v0 = cvx(randi(numel(cvx)));  d = distances(rheome.geom.edgegraph(Sx), v0)';   % plant on cortex
        f = exp(-d.^2 / (2 * (opts.SigmaMM/1e3)^2));
        e = cross(n(v0,:), [0 0 1]);  if norm(e) < 0.3, e = cross(n(v0,:), [1 0 0]); end
        e = e / norm(e);  t = e - (n * e') .* n;                       % tangential, |t| <= 1, smooth
        lp = @(x) x;
        if ~isempty(s.lbo)                                             % keep wavelengths > LowpassMM
            k = s.lbo.Lambda(:) <= (2*pi / (opts.LowpassMM/1e3))^2;  P = s.lbo.Phi(:, k);
            lp = @(x) P * (P' * (s.lbo.Mass * x));
        end
        nrm = @(x) sqrt(a' * x.^2);  emb = @(w) reshape(w.', [], 1);
        for th = opts.TiltsDeg(:)'
            Jn = f * cosd(th);  J3 = Jn .* n + (f * sind(th)) .* t;    % J . n = Jn exactly (t is tangential)
            Dam = rheome.differential.divergence(emb(J3), Sx, fg);
            Dtn = rheome.differential.divergence(emb(J3 - Jn .* n), Sx, fg);
            C = Dam - Dtn;  K = 2 * Hk .* Jn;
            cn = rheome.differential.curl(emb(Jn .* n), Sx, fg);
            r = @(x, y) i_r(x, y, a);
            if isempty(s.lbo), e60 = NaN;  r60 = NaN;
            else, e60 = nrm(lp(C) - lp(K)) / nrm(lp(K));  r60 = r(lp(C), lp(K)); end
            row = {s.surface, s.hemi, th, opts.SigmaMM, v0, nV, Hrel, nrm(C - K) / nrm(K), r(C, K), e60, r60, ...
                   nrm(C) / nrm(Dam), nrm(cn) / nrm(C)};
            X = [X; row]; %#ok<AGROW>
            b = s.surface + iif(s.hemi == "", "", "-" + s.hemi) + "/" + th + " deg";
            if cosd(th) > 1e-9
                m = ["couple_relerr" "couple_r" "normal_share"];  v = [row{8} row{9} row{12}];
                if ~isempty(s.lbo), m = [m "couple_relerr_60mm"]; v = [v e60]; end %#ok<AGROW>
                T = [T; rheome.scale.rows("ambient", m, v, "ratio", b)]; %#ok<AGROW>
            end
            if th == opts.ExampleTiltDeg
                Ex.(matlab.lang.makeValidName(char(s.surface + s.hemi))) = struct('surface', char(s.surface), ...
                    'hemi', char(s.hemi), 'tilt_deg', th, 'vertices_mm', 1e3*V, 'faces', F, 'normals', n, ...
                    'H_per_m', Hk, 'f', f, 'J', J3, 'D_amb', Dam, 'D_tan', Dtn, 'C', C, 'K', K, ...
                    'units', 'J arbitrary (peak 1); D, C, K in 1/m x J');
            end
        end
    end
    X.Properties.VariableNames = {'surface','hemi','tilt_deg','sigma_mm','centre_vertex','n_vertices','H_relerr', ...
        'couple_relerr','couple_r','couple_relerr_60mm','couple_r_60mm','normal_share','curl_normal_rel'};
end

function r = i_r(x, y, a)
% area-weighted Pearson r
    x = x - (a' * x) / sum(a);  y = y - (a' * y) / sum(a);
    r = (a' * (x .* y)) / sqrt((a' * x.^2) * (a' * y.^2));
end

function z = iif(c, a, b)
    if c, z = a; else, z = b; end
end

% Author: Diellor Basha, 2026
