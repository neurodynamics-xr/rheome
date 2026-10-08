classdef tInverseResolution < matlab.unittest.TestCase
% Resolution of an inverse, reduced to lengths. Built on a coarse sphere so the whole class
% runs in seconds: the geometry is not the point, the identities are.
%
% Author: Diellor Basha, 2026

    properties
        S; lbo; Gain; K; nV
    end

    methods (TestClassSetup)
        function build(t)
            [V, F] = i_icosphere(4);                      % 2562 vertices, h ~ 4.9 mm
            V = 0.07 * V;                                 % 70 mm radius, head-sized
            [L, M] = rheome.operators.laplace_beltrami(V, F);
            t.S = struct('Vertices', V, 'Faces', F, 'nV', size(V,1));
            t.lbo = rheome.eigen.modes(L, M, 120);   % .Phi .Lambda .Mass -- what resolution wants
            t.nV = size(V,1);
            % a crude forward model: 60 "sensors" on a shell, each reading the radial
            % component with a 1/r^2 falloff. Enough to be rank-deficient and depth-biased.
            rng(3);
            d = randn(60,3);  d = 0.11 * d ./ vecnorm(d,2,2);
            G = zeros(60, 3*t.nV);
            for c = 1:60
                dv = V - d(c,:);  r = vecnorm(dv, 2, 2);
                wgt = 1 ./ max(r, 1e-3).^2;
                G(c, :) = reshape((wgt .* (dv ./ r))', 1, []);
            end
            t.Gain = G;
            R = rheome.inverse.mne(G, struct('NoiseCov', eye(60)), ...
                    struct('ChannelTypes', {repmat({'MEG'}, 1, 60)}, ...
                           'InverseMeasure', 'amplitude', 'nVert', t.nV, 'SnrFixed', 3));
            t.K = R.ImagingKernel;
        end
    end

    methods (Test)
        function theTraceOfTheResolutionMatrixIsTheSumOfTheWienerGains(t)
            % ⭐ the identity the whole design rests on: R = VL*diag(g.*s)*VL', so trace(R)
            % is the sum of the Wiener gains -- the effective rank the regularisation sets.
            R = rheome.inverse.mne(t.Gain, struct('NoiseCov', eye(60)), ...
                    struct('ChannelTypes', {repmat({'MEG'}, 1, 60)}, ...
                           'InverseMeasure', 'amplitude', 'nVert', t.nV, 'SnrFixed', 3));
            Lam = R.SNR / mean(R.SL.^2);
            gs  = (Lam*R.SL.^2) ./ (Lam*R.SL.^2 + 1);
            % trace(K*Gain) without forming the [3nV x 3nV] product
            tr = sum(sum(t.K .* R.ImagingKernel*0, 2)) + sum(dot(t.K, t.Gain', 2));
            t.verifyEqual(tr, sum(gs), 'RelTol', 1e-8);
        end

        function moreRegularisationMeansFewerDegreesOfFreedom(t)
            g = @(snr) i_dof(t.Gain, t.nV, snr);
            t.verifyLessThan(g(1), g(3));
            t.verifyLessThan(g(3), g(30));
        end

        function theSpreadStatisticsAreOrderedAndPositive(t)
            out = rheome.inverse.resolution(t.K, t.Gain, t.lbo, Vertices=t.S.Vertices, ...
                      Faces=t.S.Faces, NumSeeds=25, Modes=0);
            t.verifyTrue(all(out.r50 > 0));
            t.verifyTrue(all(out.sd >= out.r50));         % sqrt(E[d^2]) >= the median radius
            t.verifyTrue(all(out.ple >= 0));
            t.verifyEqual(out.offHemi, zeros(25,1), 'AbsTol', 1e-12);   % one component: no leak
        end

        function theWavelengthConversionIsLinearThroughTheOrigin(t)
            % ⭐ r50 = c*wavelength is MEASURED, not assumed, and it had better be a line.
            out = rheome.inverse.resolution(t.K, t.Gain, t.lbo, Vertices=t.S.Vertices, ...
                      Faces=t.S.Faces, NumSeeds=12, Modes=0, ...
                      ... % ⚠ above the BASIS floor (40 mm here) and under the geodesic
                      ... % diameter (232 mm): outside that window r50 is not linear in it
                      Apertures=[80 110 150 200]*1e-3);
            t.verifyGreaterThan(out.r50R2, 0.99);
            t.verifyGreaterThan(out.r50PerWavelength, 0.1);
            t.verifyLessThan(out.r50PerWavelength, 0.4);
            t.verifyEqual(out.lambdaLoc, median(out.r50)/out.r50PerWavelength, 'RelTol', 1e-12);
        end

        function aWiderApertureGivesAProportionallyWiderR50(t)
            out = rheome.inverse.resolution(t.K, t.Gain, t.lbo, Vertices=t.S.Vertices, ...
                      Faces=t.S.Faces, NumSeeds=8, Modes=0, Apertures=[80 160]*1e-3);
            r = out.calibration.r50;
            t.verifyEqual(r(2)/r(1), 2, 'RelTol', 0.15);
        end

        function theMtfSumsToTheSubspaceRankAndIsNonNegative(t)
            out = rheome.inverse.resolution(t.K, t.Gain, t.lbo, Vertices=t.S.Vertices, ...
                      Faces=t.S.Faces, Normals=i_vnormals(t.S), NumSeeds=0, Modes=1:60);
            % R is PSD, so every diagonal entry in an orthonormal family is >= 0 up to the
            % round-off of an SVD-built kernel -- tolerance relative to the curve's own scale
            t.verifyGreaterThanOrEqual(min(out.mtf), -1e-5 * max(out.mtf));
            t.verifyEqual(out.dofSubspace, sum(out.mtf), 'RelTol', 1e-12);
            t.verifyLessThan(out.dofSubspace, i_dof(t.Gain, t.nV, 3) + 1e-6);   % a subspace of it
        end

        function theConstantModeIsKeptInTheRankButNotInTheBands(t)
            % ⚠ the regression: 2*pi/sqrt(1e-12) is 6e6 m and it used to set the top octave.
            out = rheome.inverse.resolution(t.K, t.Gain, t.lbo, Vertices=t.S.Vertices, ...
                      Faces=t.S.Faces, Normals=i_vnormals(t.S), NumSeeds=0, Modes=1:60);
            t.verifyEqual(out.mtfBandExcluded, 1);
            t.verifyEqual(numel(out.mtf), 60);            % still counted in .mtf
            t.verifyLessThan(out.mtfBandEdges(1), 1);     % top edge under a metre
            t.verifyTrue(all(isfinite(out.mtfBandEdges)));
        end

        function aMismatchedKernelAndLeadfieldAreRefused(t)
            t.verifyError(@() rheome.inverse.resolution(t.K, t.Gain(1:30,:), t.lbo, ...
                Vertices=t.S.Vertices, Faces=t.S.Faces, NumSeeds=2, Modes=0), ...
                'inverse:resolution:chan');
            t.verifyError(@() rheome.inverse.resolution(t.K, t.Gain, t.lbo, ...
                Vertices=t.S.Vertices, Faces=t.S.Faces, GlobalIdx=1:5, NumSeeds=2, Modes=0), ...
                'inverse:resolution:idx');
        end
    end
end

function d = i_dof(G, nV, snr)
    R = rheome.inverse.mne(G, struct('NoiseCov', eye(size(G,1))), ...
            struct('ChannelTypes', {repmat({'MEG'}, 1, size(G,1))}, ...
                   'InverseMeasure', 'amplitude', 'nVert', nV, 'SnrFixed', snr));
    L = R.SNR / mean(R.SL.^2);
    d = sum((L*R.SL.^2) ./ (L*R.SL.^2 + 1));
end

function N = i_vnormals(S)
    P = S.Vertices;  F = S.Faces;
    fn = cross(P(F(:,2),:) - P(F(:,1),:), P(F(:,3),:) - P(F(:,1),:), 2);
    N = zeros(size(P,1), 3);
    for c = 1:3, N(:,c) = accumarray(F(:), repmat(fn(:,c), 3, 1), [size(P,1) 1]); end
    N = N ./ max(vecnorm(N, 2, 2), eps);
end

function [V, F] = i_icosphere(n)
    p = (1+sqrt(5))/2;
    V = [-1 p 0; 1 p 0; -1 -p 0; 1 -p 0; 0 -1 p; 0 1 p; 0 -1 -p; 0 1 -p; p 0 -1; p 0 1; -p 0 -1; -p 0 1];
    V = V ./ vecnorm(V, 2, 2);
    F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12; 2 6 10; 6 12 5; 12 11 3; 11 8 7; 8 2 9; ...
         4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10; 5 10 6; 3 5 12; 7 3 11; 9 7 8; 10 9 2];
    for k = 1:n
        nF = [];  M = containers.Map('KeyType','char','ValueType','double');
        for i = 1:size(F,1)
            a = F(i,1); b = F(i,2); c = F(i,3);
            [V, M, ab] = i_mid(V, M, a, b);  [V, M, bc] = i_mid(V, M, b, c);
            [V, M, ca] = i_mid(V, M, c, a);
            nF = [nF; a ab ca; b bc ab; c ca bc; ab bc ca];   %#ok<AGROW>
        end
        F = nF;
    end
end

function [V, M, m] = i_mid(V, M, a, b)
    k = sprintf('%d_%d', min(a,b), max(a,b));
    if M.isKey(k), m = M(k); return; end
    v = (V(a,:) + V(b,:)) / 2;  v = v / norm(v);
    V = [V; v];  m = size(V,1);  M(k) = m;
end
% Author: Diellor Basha, 2026
