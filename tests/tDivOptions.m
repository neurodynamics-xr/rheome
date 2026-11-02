classdef tDivOptions < matlab.unittest.TestCase
    % rheome.differential divergence / curl options (tan, 3d, curv), the mean curvature and 2 H J_n against the goldens nxr wrote
    % (tests/goldens/div-options, synthetic: ellipsoid 70/55/45 mm, icosphere 4; PR neurodynamics-xr/nxr-analysis#11), and the
    % curvature term through the Dirac.
    %
    % Author: Diellor Basha, 2026

    properties
        G      % struct of golden arrays
        S      % surface in metres, positions cast to single as nxr does
        fg
        J      % [3nV x 1]
        R      % divcurl result
    end

    methods (TestClassSetup)
        function load(t)
            d = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'div-options');
            M = jsondecode(fileread(fullfile(d, 'manifest.json')));
            nm = fieldnames(M.arrays);
            for i = 1:numel(nm)
                a = M.arrays.(nm{i}); sh = double(a.shape(:)).';
                f = fopen(fullfile(d, a.file), 'r'); x = fread(f, prod(sh), ['*' strrep(a.dtype, 'float64', 'double')]); fclose(f);
                if numel(sh) == 1, x = x(:); else, x = permute(reshape(x, fliplr(sh)), numel(sh):-1:1); end
                t.G.(nm{i}) = double(x);
            end
            t.S = struct('Vertices', double(single(t.G.V_mm / 1000)), 'Faces', t.G.F + 1);   % nxr runs on float32 metres
            t.fg = rheome.operators.face_gradient(t.S.Vertices, t.S.Faces);
            t.J = reshape(t.G.J.', [], 1);
            t.R = rheome.differential.divcurl(t.J, t.S, t.fg);
        end
    end

    methods (Test)
        function everyGoldenMap(t)
            for k = ["div_tan" "div_curv" "div3d" "curl_tan" "curl_curv" "curl3d"]
                ref = t.G.(k);  e = max(abs(t.R.(k) - ref)) / max(abs(ref));
                t.verifyLessThan(e, 1e-6, char(k));
                fprintf('  %-10s max rel %.2e\n', k, e);
            end
        end

        function meanCurvatureAndTwoHJn(t)
            H = rheome.differential.meancurvature(t.S, t.fg);
            e = max(abs(H - t.G.H)) / max(abs(t.G.H));  fprintf('  H          max rel %.2e\n', e);
            t.verifyLessThan(e, 1e-6);
            c = rheome.differential.curvaturesource(t.J, t.S, t.fg);
            e = max(abs(c - t.G.div_curv_2H)) / max(abs(t.G.div_curv_2H));  fprintf('  2 H J_n    max rel %.2e\n', e);
            t.verifyLessThan(e, 1e-6);
        end

        function unprojectedPositionsAlsoMatch(t)
            % Rheome does not cast positions to single (nxr and the Editor do): the effect on the maps
            S0 = t.S;  S0.Vertices = t.G.V_mm / 1000;  R0 = rheome.differential.divcurl(t.J, S0);
            for k = ["div_tan" "div_curv" "div3d" "curl_tan" "curl_curv" "curl3d"]
                e = max(abs(R0.(k) - t.G.(k))) / max(abs(t.G.(k)));  fprintf('  %-10s double positions, max rel %.2e\n', k, e);
                t.verifyLessThan(e, 2e-5, char(k));
            end
        end

        function identityToRounding(t)
            t.verifyLessThan(max(abs(t.R.div3d - t.R.div_tan - t.R.div_curv)) / max(abs(t.R.div3d)), 1e-12);
            t.verifyLessThan(max(abs(t.R.curl3d - t.R.curl_tan - t.R.curl_curv)) / max(abs(t.R.curl3d)), 1e-12);
            % a rough mesh, random field
            rng(0);  [V, F] = rheome.geom.icosphere(3);  V = 0.06 * V .* [1.2 1 0.8] + 1e-3 * randn(size(V));
            S = struct('Vertices', V, 'Faces', F);  Jr = randn(3*size(V,1), 3);  M = rheome.differential.divcurl(Jr, S);
            t.verifyLessThan(max(abs(M.div3d - M.div_tan - M.div_curv), [], 'all') / max(abs(M.div3d), [], 'all'), 1e-12);
            t.verifyLessThan(max(abs(M.curl3d - M.curl_tan - M.curl_curv), [], 'all') / max(abs(M.curl3d), [], 'all'), 1e-12);
        end

        function defaultsUnchanged(t)
            f = t.fg;  J = t.J;
            t.verifyEqual(rheome.differential.curl(J, t.S, f), t.R.curl_tan);
            t.verifyEqual(rheome.differential.curl(J, t.S, f, true), t.R.curl_tan);
            t.verifyEqual(rheome.differential.curl(J, t.S, f, false), t.R.curl3d);
            t.verifyEqual(rheome.differential.divergence(J, t.S, f), t.R.div3d);
            t.verifyEqual(rheome.differential.divergence(J, t.S, f, false), t.R.div3d);
            t.verifyEqual(rheome.differential.divergence(J, t.S, f, true), t.R.div_tan);
            t.verifyEqual(rheome.differential.divergence(J, t.S, f, "tan"), t.R.div_tan);
            t.verifyError(@() rheome.differential.divergence(J, t.S, f, "x"), 'differential:divergence:option');
            t.verifyError(@() rheome.differential.curl(J, t.S, f, "x"), 'differential:curl:option');
        end

        function sphereNormalCurrent(t)
            % J = c n on a sphere of radius R: div_tan = 0, div3d = div_curv = 2c/R, curl_curv ~ 0, H = 1/R; converges with refinement
            Rr = 0.06;  c0 = 2;  e = zeros(3, 4);
            for l = 3:5
                [V, F] = rheome.geom.icosphere(l);  S = struct('Vertices', Rr*V, 'Faces', F);  fg = rheome.operators.face_gradient(S.Vertices, F);
                n = fg.W * fg.FaceNormal;  n = n ./ vecnorm(n, 2, 2);  J = reshape((c0*n).', [], 1);  tgt = 2*c0/Rr;
                M = rheome.differential.divcurl(J, S, fg);  H = rheome.differential.meancurvature(S, fg);
                rms = @(x, r) sqrt(mean((x - r).^2)) / tgt;
                e(l-2, :) = [max(abs(M.div_tan))/tgt, rms(M.div3d, tgt), rms(rheome.differential.curvaturesource(J, S, fg), tgt), ...
                             sqrt(mean(M.curl_curv.^2)) / (c0/Rr)];
                t.verifyLessThan(sqrt(mean((H - 1/Rr).^2)) * Rr, 0.1);
            end
            fprintf('  sphere ico3-5: max|div_tan|/t, rms div3d, rms 2HJn, rms curl_curv (relative)\n'); disp(e)
            t.verifyLessThan(e(:,1), 1e-12);
            t.verifyLessThan(e(3,2:4), e(1,2:4));
            t.verifyLessThan(e(3,2), 0.025);  t.verifyLessThan(e(3,3), 0.02);  t.verifyLessThan(e(3,4), 0.02);
        end

        function diracReadsTheCurvatureTerm(t)
            % (1) the intrinsic Dirac of the NORMAL part reads div_curv / curl_curv exactly (it reads the ambient operators)
            n = t.fg.W * t.fg.FaceNormal;  n = n ./ vecnorm(n, 2, 2);  Jv = reshape(t.J, 3, []).';
            Jn = sum(Jv .* n, 2);  Jc = reshape((Jn .* n).', [], 1);
            P = rheome.differential.diracparts(Jc, t.S);
            t.verifyLessThan(max(abs(P.div - t.G.div_curv)) / max(abs(t.G.div_curv)), 1e-6);
            t.verifyLessThan(max(abs(P.curl - t.G.curl_curv)) / max(abs(t.G.curl_curv)), 1e-6);
            P3 = rheome.differential.diracparts(t.J, t.S);
            t.verifyLessThan(max(abs(P3.div - t.G.div3d)) / max(abs(t.G.div3d)), 1e-6);   % so the Dirac div is div_tan + div_curv
            % (2) extrinsic - intrinsic: D_HY - D_ext = H_f (integrated mean curvature per face, dihedral angles, an estimator
            %     independent of the cotan Laplacian). On q = J_n n:  div_curv_f ~ 2 (H_f / A_f) (q . n_f) = 2 H J_n
            nF = size(t.S.Faces, 1);  F = t.S.Faces;  D = rheome.operators.spin_dirac(t.S.Vertices, F);
            dif = D.intrinsicHY - D.extrinsic;
            t.verifyEqual(nnz(dif - spdiags(diag(dif), 0, 4*nF, 4*nF)), 0);                    % a pure multiplication by H_f
            qv = reshape(Jc, 3, []).';  qf = (qv(F(:,1),:) + qv(F(:,2),:) + qv(F(:,3),:)) / 3;
            q = zeros(4*nF, 1);  q(2:4:end) = qf(:,1);  q(3:4:end) = qf(:,2);  q(4:4:end) = qf(:,3);
            Y = dif * q;  yd = sum([Y(2:4:end) Y(3:4:end) Y(4:4:end)] .* t.fg.FaceNormal, 2) ./ t.fg.FaceArea;
            tof = @(x) (x(F(:,1)) + x(F(:,2)) + x(F(:,3))) / 3;
            ref = tof(t.G.div_curv);  e2 = norm(2*yd - ref) / norm(ref);  e1 = norm(yd - ref) / norm(ref);
            fprintf('  Dirac route: 2 (D_HY - D_ext) q.n/A vs face-mean div_curv rel L2 %.3e (factor 1: %.3e)\n', e2, e1);
            t.verifyLessThan(e2, 0.05);  t.verifyGreaterThan(e1, 0.4);
            c2 = tof(t.G.div_curv_2H);  fprintf('  face-mean div_curv vs 2HJn rel L2 %.3e\n', norm(ref - c2) / norm(c2));
        end
    end
end
