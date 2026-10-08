classdef tDiracParts < matlab.unittest.TestCase
% The intrinsic Dirac applied to a vector field IS its curl, divergence and normal-component gradient,
% exactly; and on a tau = 0 Dirac basis that split holds band by band.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function theThreePartsAreCurlDivAndNormalGradient(tc)
            [V, F] = rheome.geom.icosphere(3);  rng(1);
            V = V .* (1 + 0.2*sin(3*V(:,1)).*cos(2*V(:,2)));          % not a sphere
            for FF = {F, F(:, [1 3 2])}                                % both orientations
                S = struct('Vertices', V, 'Faces', FF{1});
                J = randn(3*size(V,1), 2);
                P = rheome.differential.diracparts(J, S);
                tc.verifyEqual(P.curl, rheome.differential.curl(J, S), 'AbsTol', 1e-12 * max(abs(P.curl), [], 'all'));
                tc.verifyEqual(P.div, rheome.differential.divergence(J, S), 'AbsTol', 1e-12 * max(abs(P.div), [], 'all'));
                G = i_gradnormal(J(:,1), V, FF{1});
                tc.verifyEqual(reshape(P.gradNF(:,1), 3, [])', G, 'AbsTol', 1e-12 * max(abs(G), [], 'all'));
            end
        end

        function theSplitIsTheDiracEnergy(tc)
            [V, F] = rheome.geom.icosphere(3);  rng(2);  S = struct('Vertices', V, 'Faces', F);
            J = randn(3*size(V,1), 3);  P = rheome.differential.diracparts(J, S);
            L = rheome.operators.dirac_intrinsic_sq(V, F);  Q = i_quat(J);
            tc.verifyEqual(P.energy.total, diag(Q' * L * Q)', 'RelTol', 1e-12);
        end

        function aGradientIsDivergentAndARotatedGradientIsRotational(tc)
            [V, F] = rheome.geom.icosphere(4);  S = struct('Vertices', V, 'Faces', F);  n = V ./ vecnorm(V, 2, 2);
            phi = V(:,3) + 0.5*V(:,1).*V(:,2);
            g = [0.5*V(:,2), 0.5*V(:,1), ones(size(V,1),1)];                 % ambient grad of phi
            gt = g - sum(g.*n, 2) .* n;                                        % its tangential part
            Pg = rheome.differential.diracparts(reshape(gt', [], 1), S);
            Pr = rheome.differential.diracparts(reshape(cross(n, gt, 2)', [], 1), S);
            Pn = rheome.differential.diracparts(reshape((phi .* n)', [], 1), S);
            tc.verifyLessThan(Pg.energy.curl / Pg.energy.total, 1e-4);        % a gradient has no curl
            tc.verifyLessThan(Pr.energy.div / Pr.energy.total, 1e-4);         % a rotated gradient no div
            tc.verifyLessThan(Pn.energy.curl / Pn.energy.total, 1e-4);        % a normal field no curl
            % ⚠ a TANGENT field still has a third part: grad(X . n_f) with n_f frozen per face is the shape
            % operator acting on X, so on the unit sphere (S = I) its energy is the integral of |X|^2
            M = rheome.operators.mass(V, F, 'galerkin');
            tc.verifyEqual(Pg.energy.gradN, sum(M * sum(gt.^2, 2)), 'RelTol', 5e-3);
            tc.verifyEqual(Pr.energy.gradN, Pg.energy.gradN, 'RelTol', 5e-3);   % both -> int |X|^2
        end

        function onATauZeroBasisTheSplitHoldsBandByBand(tc)
            [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', V, 'Faces', F);  rng(3);
            db = rheome.eigen.dirac_frame(V, F, 0, 160, [], {(1:size(V,1))'}, false);   % tau 0, physical lambda
            Phi = full(db.Phi);  lam = db.Lambda;
            c = randn(numel(lam), 1);  x = Phi * c;
            P = rheome.differential.diracparts(x, S);
            tc.verifyEqual(P.energy.total, lam' * c.^2, 'RelTol', 1e-8);          % ||D x||^2 = sum lambda c^2
            gfb = rheome.graphfilterbank(lam(lam > 1e-9*max(lam)), 'Wavelet', 'logitersine', 'VoicesPerOctave', 1);
            G = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:gfb.NumMembers, 'uni', 0));
            tc.verifyEqual(sum(G.^2, 2), ones(size(lam)), 'AbsTol', 1e-6);         % tight on this spectrum
            Eb = zeros(gfb.NumMembers, 3);
            for m = 1:gfb.NumMembers
                Pm = rheome.differential.diracparts(Phi * (G(:,m) .* c), S);
                Eb(m, :) = [Pm.energy.curl Pm.energy.div Pm.energy.gradN];
                tc.verifyEqual(sum(Eb(m,:)), lam' * (G(:,m).^2 .* c.^2), 'RelTol', 1e-8);
            end
            tc.verifyEqual(sum(Eb, 'all'), P.energy.total, 'RelTol', 1e-6);
        end
    end
end

function Q = i_quat(J)
    nV = size(J,1)/3;  Q = zeros(4*nV, size(J,2));  Q(setdiff(1:4*nV, 1:4:4*nV), :) = J;
end

% grad (X . n_f) per face, from face_gradient's hat gradients and the face normal
function G = i_gradnormal(J, V, F)
    fg = rheome.operators.face_gradient(V, F);  X = reshape(J, 3, [])';  N = fg.FaceNormal;  nF = size(F,1);
    G = zeros(nF, 3);
    for c = 1:3
        i = F(:,c);  k = sub2ind(size(fg.Gx), (1:nF)', i);
        G = G + [full(fg.Gx(k)) full(fg.Gy(k)) full(fg.Gz(k))] .* sum(X(i,:) .* N, 2);
    end
end

% Author: Diellor Basha, 2026
