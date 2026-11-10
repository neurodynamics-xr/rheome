classdef tRegister < matlab.unittest.TestCase
% rheome.register (hemispheres, registration, shepard_weights, mode_registration, check_mode_registration, to_glued / to_input,
% register_values, nearest_vertices) and rheome.inverse.hodge_register on a synthetic
% two-piece surface (two icospheres side by side, a medial cap cut on each, each its own registered sphere), registered onto
% itself: what must hold exactly. The sampling is held to nxr's own (nxr.register.sphere) on a synthetic sphere when nxr and uv are
% found ($NXR, the nxr checkout; skipped otherwise). The parity with the editor's stored registration is in tNxrStoreParity.
%
% Author: Diellor Basha, 2026

    properties
        side, reg
    end

    methods (TestClassSetup)
        function build(tc)
            [V0, F0] = rheome.geom.icosphere(3); n0 = size(V0, 1);
            P = [V0 - [1.5 0 0]; V0 + [1.5 0 0]]; F = [F0; F0 + n0];
            wall = [V0(:, 1) > 0.75; -V0(:, 1) > 0.75];
            [X, G, keep] = rheome.geom.glue(P, F, wall);
            H = rheome.operators.hodge_dirac(X, G);
            pkg = rheome.eigen.hodge_package(rheome.eigen.hodge_dirac(H, 40), H);
            tc.side = rheome.register.solve(pkg, keep, X);
            tc.reg = rheome.register.surfaces(P, F, [V0; V0], P, F, [V0; V0], SubjectLeft="-x", TemplateLeft="-x");
        end
    end

    methods (Test)
        function hemispheresByTheAxis(tc)
            P = tc.reg.subject_positions; F = tc.reg.subject_faces;
            h = rheome.register.hemispheres(P, F, "-x");
            tc.verifyTrue(all(P(h == 'L', 1) < 0) && all(P(h == 'R', 1) > 0) && ~any(h == '?'));
            tc.verifyEqual(rheome.register.hemispheres(P, F, "+x"), char('L' + 'R' - h), 'the other axis swaps them');
        end

        function registrationIsAPartitionOfUnity(tc)
            [T, R] = rheome.register.mode_registration(tc.side, tc.side, tc.reg);
            m = tc.side.kept;
            tc.verifyEqual(full(sum(R, 2)), ones(size(R, 1), 1), 'AbsTol', 1e-12, 'every template vertex a blend of subject vertices');
            tc.verifyLessThan(full(max(abs(R(1:m, 1:m) - speye(m)), [], 'all')), 1e-12, 'cortex onto itself: the identity');
            tc.verifyEqual(full(max(abs(R(:, m+1:end)), [], 'all')), 0, 'the subject''s bridge is never sampled');
            k = tc.side.kind(:);
            tc.verifyEqual(max(abs(T(k == 0, k == 2)), [], 'all') + max(abs(T(k == 2, k == 0)), [], 'all'), 0, 'block-diagonal by kind');
        end

        function registeringWeightsRegistersTheMaps(tc)
            [T, R] = rheome.register.mode_registration(tc.side, tc.side, tc.reg);
            pkg = tc.side; K = numel(pkg.mu);
            B = randn(RandStream('mt19937ar', 'Seed', 5), 2 * K, 6); x = randn(6, 3);
            Bt = rheome.inverse.hodge_register(B, T, pkg.mu, pkg.mu); w = B * x;
            tc.verifyEqual(Bt * x, [T * w(1:K, :); (T * (pkg.mu .* w(K+1:end, :))) ./ pkg.mu], 'AbsTol', 1e-10);
            tc.verifyEqual(rheome.inverse.hodge_register(B, T), [T * B(1:K, :); T * B(K+1:end, :)], 'AbsTol', 1e-12);
            tc.verifyError(@() rheome.inverse.hodge_register(B, T, pkg.mu), 'rheome:hodge_register');
            C = rheome.register.check_mode_registration(T, R, tc.side, tc.side, tc.reg);
            tc.verifyEqual(height(C), 15); tc.verifyTrue(all(C.l == repelem((1:3)', 2 * (1:3)' + 1)));
        end

        function shepardWeightsAreBrainstorms(tc)
            % k nearest, d2 (0 -> eps), R2 the k-th, ((R2 - d2) / (R2 d2))^p over the first k - 1, normalised; a point on a vertex: weight 1
            X = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
            Y = [1 0 0; 0.6 0.8 0];
            [idx, w] = rheome.register.shepard_weights(X, Y, 3, 2);
            tc.verifyEqual(size(idx), [2 2]); tc.verifyEqual(idx(1, 1), 1); tc.verifyEqual(w(1, :), [1 0], 'AbsTol', 1e-12, 'overlapping a vertex: its own');
            d2all = sort(sum((X - Y(2, :)) .^ 2, 2)); d2 = d2all(1:2)'; R2 = d2all(3);
            tc.verifyEqual(idx(2, :), [2 1]);
            e = ((R2 - d2) ./ (R2 * d2)) .^ 2; tc.verifyEqual(w(2, :), e / sum(e), 'AbsTol', 1e-12);
            [i1, w1] = rheome.register.shepard_weights(X, Y, 1); tc.verifyEqual([i1(:, 1) w1], [1 1; 2 1]);
        end

        function shepardRegistrationIsNormalisedOnTheHemisphere(tc)
            [V0, F0] = rheome.geom.icosphere(2); Y = randn(RandStream('mt19937ar', 'Seed', 3), 30, 3);
            h = repmat('L', size(V0, 1), 1);
            R = rheome.register.registration(V0, F0, h, Y, repmat('L', 30, 1), Sampling="shepard", Neighbours=5, Power=2);
            tc.verifyEqual(full(sum(R, 2)), ones(30, 1), 'AbsTol', 1e-12); tc.verifyEqual(full(max(sum(R ~= 0, 2))), 4);
            tc.verifyError(@() rheome.register.registration(V0, F0, h, Y, h(1:30), Sampling="spline"), 'MATLAB:validators:mustBeMember');
        end

        function valuesAndPointsByVertices(tc)
            [V0, F0] = rheome.geom.icosphere(2); n = size(V0, 1); h = repmat('L', n, 1);
            R = rheome.register.registration(V0, F0, h, V0, h, Sampling="nearest");
            x = (1:n)'; tc.verifyEqual(rheome.register.register_values(R, [x NaN(n, 1)]), [x zeros(n, 1)]);
            R(5, :) = 0; y = rheome.register.register_values(R, x); tc.verifyTrue(isnan(y(5)) && ~any(isnan(y([1:4 6:end]))));
            tc.verifyEqual(rheome.register.nearest_vertices(V0, F0, h, V0(10:20, :), h(10:20)), (10:20)');
        end

        function glueAndBackRoundTrips(tc)
            s = tc.side; P = tc.reg.subject_positions; x = P(:, 1) + 2 * P(:, 2);
            y = rheome.register.to_glued(s, x, P);
            tc.verifyEqual(y(1:s.kept), x(s.keep(1:s.kept)));
            tc.verifyEqual(size(y, 1), size(s.positions, 1));
            z = rheome.register.to_input(s, y, P);
            tc.verifyEqual(z(s.keep(1:s.kept)), y(1:s.kept)); tc.verifyEqual(size(z, 1), size(P, 1));
            u = rheome.register.register_in_modes(eye(numel(s.mu)), s, s, y);
            tc.verifyEqual(size(u), size(y));
        end

        function samplingIsNxrs(tc)
            nxr = getenv('NXR');
            tc.assumeTrue(~isempty(nxr) && isfile(fullfile(nxr, 'pyproject.toml')) && system('command -v uv >/dev/null') == 0, 'nxr ($NXR) or uv not found');
            [V0, F0] = rheome.geom.icosphere(2); n = size(V0, 1);
            Y = randn(RandStream('mt19937ar', 'Seed', 11), 40, 3);
            d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
            i_write(d, 'X', V0); i_write(d, 'F', F0 - 1); i_write(d, 'Y', Y);
            [st, out] = system(sprintf('uv run --frozen --project "%s" python -I "%s" "%s" 2>&1', nxr, fullfile(fileparts(mfilename('fullpath')), 'nxr_register.py'), d));
            tc.assertEqual(st, 0, out);
            h = repmat('L', n, 1); hd = repmat('L', 40, 1);
            for s = ["barycentric" "nearest" "shepard"]
                R = rheome.register.registration(V0, F0, h, Y, hd, Sampling=s, Neighbours=6, Power=3);
                tc.verifyEqual(full(R), reshape(i_read(d, "R_" + s), n, [])', 'AbsTol', 1e-12, char(s));
            end
            [idx, w] = rheome.register.shepard_weights(V0 ./ vecnorm(V0, 2, 2), Y ./ vecnorm(Y, 2, 2), 6, 3);
            tc.verifyEqual(idx, reshape(i_read(d, "S_idx"), 5, [])' + 1); tc.verifyEqual(w, reshape(i_read(d, "S_w"), 5, [])', 'AbsTol', 1e-12);
        end
    end
end

function i_write(d, name, x)
fid = fopen(fullfile(d, name + ".f64"), 'w', 'ieee-le'); fwrite(fid, x', 'double'); fclose(fid);
end

function x = i_read(d, name)
fid = fopen(fullfile(d, name + ".f64"), 'r', 'ieee-le'); x = fread(fid, Inf, 'double'); fclose(fid);
end
