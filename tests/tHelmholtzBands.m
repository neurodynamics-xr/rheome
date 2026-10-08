classdef tHelmholtzBands < matlab.unittest.TestCase
% Source and vortex content of a field, band by band, through its Helmholtz potentials: the band
% energies are exact, a source lands in the irrotational part and a vortex in the solenoidal part, and a
% potential planted in one band is found in that band.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function bandEnergiesAddUpToTheDirichletEnergy(tc)
            [S, lbo] = i_sphere();  X = i_atoms(S, lbo, 2);
            B = rheome.differential.helmholtzbands(X, S, lbo);
            tc.verifyEqual(B.capture, ones(2, 2), 'AbsTol', 1e-9);              % complete basis
            tc.verifyEqual(sum(B.Eirr, 1), B.irrTotal, 'RelTol', 1e-9);
            tc.verifyEqual(sum(B.Esol, 1), B.solTotal, 'RelTol', 1e-9);
        end

        function aSourceIsIrrotationalAndAVortexSolenoidal(tc)
            [S, lbo] = i_sphere();  X = i_atoms(S, lbo, 2);
            B = rheome.differential.helmholtzbands(X, S, lbo);
            tc.verifyLessThan(sum(B.Esol(:,1)) / sum(B.Eirr(:,1)), 0.01);       % source
            tc.verifyLessThan(sum(B.Eirr(:,2)) / sum(B.Esol(:,2)), 0.01);       % vortex
        end

        function aPotentialPlantedInOneBandIsFoundThere(tc)
            [S, lbo] = i_sphere();
            gfb = rheome.graphfilterbank(lbo.Lambda, 'Wavelet', 'logitersine', 'VoicesPerOctave', 2);
            for m0 = [4 6 8]
                X = i_atoms(S, lbo, 2, feval(gain(gfb, m0), lbo.Lambda));
                B = rheome.differential.helmholtzbands(X, S, lbo);
                [~, mi] = max(B.Eirr(:,1));  [~, ms] = max(B.Esol(:,2));
                tc.verifyLessThanOrEqual(abs([mi ms] - m0), [1 1], sprintf('planted band %d', m0));
            end
        end
    end
end

% icosphere(3) with its COMPLETE Laplace-Beltrami eigenbasis, so the bank is exactly tight
function [S, lbo] = i_sphere()
    [V, F] = rheome.geom.icosphere(3);
    S = struct('Vertices', V, 'Faces', F, 'VertNormals', V ./ vecnorm(V, 2, 2));
    [L, M] = rheome.operators.laplace_beltrami(V, F);
    [P, D] = eig(full(L), full(M), 'chol');  [lam, o] = sort(max(diag(D), 0));  P = P(:, o);
    P = P ./ sqrt(sum(P .* (M * P), 1));
    lbo = struct('Phi', P, 'Lambda', lam, 'Mass', M, 'L', L);
end

% [source, vortex] fields built on a potential psi at a vertex: tangential grad(psi) and n x grad(psi);
% psi = the LB kernel g(lambda) at the vertex (default a heat kernel of 0.3 radius)
function X = i_atoms(S, lbo, v0, g)
    if nargin < 4, g = exp(-lbo.Lambda * 0.3^2 / 2); end
    d = zeros(size(S.Vertices, 1), 1);  d(v0) = 1;
    psi = lbo.Phi * (g .* (lbo.Phi' * (lbo.Mass * d)));
    fg = rheome.operators.face_gradient(S.Vertices, S.Faces);  n = S.VertNormals;
    gr = [fg.W*(fg.Gx*psi) fg.W*(fg.Gy*psi) fg.W*(fg.Gz*psi)];  gr = gr - sum(gr.*n, 2) .* n;
    X = [reshape(gr', [], 1) reshape(cross(n, gr, 2)', [], 1)];
end

% Author: Diellor Basha, 2026
