classdef tFlowBandlimit < matlab.unittest.TestCase
% A current field band-limited through its Helmholtz potentials by the tight graph-wavelet bank: smooth
% where the per-vertex field is not, sign-faithful div/curl, lossless when every band is kept; and the
% Dirichlet-wavelength smoothness measure that reports it.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function aModeReadsItsOwnWavelength(tc)
            [S, lbo] = i_sphere();  k = [5 30 120];
            s = rheome.flow.smoothness(lbo.Phi(:, k), S, lbo);
            tc.verifyEqual(s.wavelengthMM, 1e3 * 2*pi ./ sqrt(lbo.Lambda(k))', 'RelTol', 1e-9);
            tc.verifyTrue(all(isnan(s.coherence)));
        end

        function parallelArrowsAreCoherent(tc)
            [S, lbo] = i_sphere();  nV = size(S.Vertices, 1);
            s = rheome.flow.smoothness(repmat([1; 0; 0], nV, 1), S, lbo);
            tc.verifyEqual(s.coherence, 1, 'AbsTol', 1e-12);
            r = randn(3*nV, 1);  s = rheome.flow.smoothness(r, S, lbo);
            tc.verifyLessThan(abs(s.coherence), 0.1);
        end

        function keepingEveryBandIsLossless(tc)
            [S, lbo] = i_sphere();  X = i_atoms(S, lbo, 2);
            B = rheome.flow.bandlimit(X, S, lbo, CutoffMM=1);
            tc.verifyEqual(B.h, ones(size(B.h)), 'AbsTol', 1e-9);
            H = rheome.differential.helmholtz(X, S);  T = H.Virr + H.Vsol;   % the tangential part, same discretisation
            tc.verifyLessThan(vecnorm(B.J - T) ./ vecnorm(T), [1e-6 1e-6]);
        end

        function noiseIsRemovedAndDivCurlKeepTheirSign(tc)
            [S, lbo] = i_sphere();  rng(1);  X = i_atoms(S, lbo, 2);
            Xn = X + 2 * std(X(:)) * randn(size(X));         % mesh-scale noise, free orientation
            B = rheome.flow.bandlimit(Xn, S, lbo, CutoffMM=700);
            raw = rheome.flow.smoothness(Xn, S, lbo);  bl = rheome.flow.smoothness(B.J, S, lbo);
            tc.verifyGreaterThan(bl.wavelengthMM, 2 * raw.wavelengthMM);
            tc.verifyGreaterThan(bl.coherence, raw.coherence + 0.2);
            dRaw = rheome.differential.divergence(Xn, S);  cRaw = rheome.differential.curl(Xn, S);
            tc.verifyGreaterThan(rheome.flow.smoothness(B.Div(:,1), S, lbo).wavelengthMM, ...
                                 2 * rheome.flow.smoothness(dRaw(:,1), S, lbo).wavelengthMM);
            % against the CLEAN field's pointwise operators: the source's div, the vortex's curl
            d0 = rheome.differential.divergence(X, S);  c0 = rheome.differential.curl(X, S);
            tc.verifyGreaterThan(corr(B.Div(:,1), d0(:,1)), 0.8);
            tc.verifyGreaterThan(corr(B.Curl(:,2), c0(:,2)), 0.8);
            tc.verifyGreaterThan(corr(B.Div(:,1), d0(:,1)), corr(dRaw(:,1), d0(:,1)));
            tc.verifyGreaterThan(corr(B.Curl(:,2), c0(:,2)), corr(cRaw(:,2), c0(:,2)));
        end

        function aCutoffAboveTheBankIsAnError(tc)
            [S, lbo] = i_sphere();  X = i_atoms(S, lbo, 2);
            tc.verifyError(@() rheome.flow.bandlimit(X, S, lbo, CutoffMM=1e9), 'flow:bandlimit:cutoff');
        end
    end
end

function [S, lbo] = i_sphere()
    [V, F] = rheome.geom.icosphere(3);
    S = struct('Vertices', V, 'Faces', F, 'VertNormals', V ./ vecnorm(V, 2, 2));
    [L, M] = rheome.operators.laplace_beltrami(V, F);
    [P, D] = eig(full(L), full(M), 'chol');  [lam, o] = sort(max(diag(D), 0));  P = P(:, o);
    P = P ./ sqrt(sum(P .* (M * P), 1));
    lbo = struct('Phi', P, 'Lambda', lam, 'Mass', M, 'L', L);
end

% [source, vortex]: tangential grad(psi) and n x grad(psi), psi a heat kernel of 0.3 radius at v0
function X = i_atoms(S, lbo, v0)
    g = exp(-lbo.Lambda * 0.3^2 / 2);
    d = zeros(size(S.Vertices, 1), 1);  d(v0) = 1;
    psi = lbo.Phi * (g .* (lbo.Phi' * (lbo.Mass * d)));
    fg = rheome.operators.face_gradient(S.Vertices, S.Faces);  n = S.VertNormals;
    gr = [fg.W*(fg.Gx*psi) fg.W*(fg.Gy*psi) fg.W*(fg.Gz*psi)];  gr = gr - sum(gr.*n, 2) .* n;
    X = [reshape(gr', [], 1) reshape(cross(n, gr, 2)', [], 1)];
end

% Author: Diellor Basha, 2026
