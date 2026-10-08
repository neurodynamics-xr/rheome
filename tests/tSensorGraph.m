classdef tSensorGraph < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (Test)

        function laplacianRowsSumToZero(tc)
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3));
            tc.verifyEqual(full(sum(G.L, 2)), zeros(36,1), 'AbsTol', 1e-12);
        end

        function laplacianIsSymmetricPositiveSemidefinite(tc)
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3));
            L = full(G.L);
            tc.verifyEqual(L, L.', 'AbsTol', 1e-14);
            tc.verifyGreaterThan(min(eig(L)), -1e-10);
        end

        function theConstantVectorIsTheZeroMode(tc)
            % A connected graph has exactly one zero eigenvalue and the constant is its
            % eigenvector. Two zeros would mean the kNN graph had fallen apart.
            G  = rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3));
            ev = sort(eig(full(G.L)));
            tc.verifyEqual(ev(1), 0, 'AbsTol', 1e-10);
            tc.verifyGreaterThan(ev(2), 1e-8);
            tc.verifyEqual(full(G.L * ones(36,1)), zeros(36,1), 'AbsTol', 1e-12);
        end

        function weightsAreGaussianInDistance(tc)
            % Two sensors one pitch apart, sigma = pitch: w = exp(-1/2). Pinning the
            % convention, because alpha in the calibration is linear in it.
            arr = rheome.sensors.grid('Size', [2 2], 'Pitch', 3e-3);
            G   = rheome.sensors.graph(arr, 'Sigma', 3e-3, 'K', 3);
            tc.verifyEqual(full(G.W(1,2)), exp(-0.5), 'RelTol', 1e-12);
            tc.verifyEqual(full(G.W(1,4)), exp(-1.0), 'RelTol', 1e-12);   % the diagonal
        end

        function theGraphIsSymmetrisedByMaxNotByAverage(tc)
            % On an 8-site chain (pitch h) with K=2: vertex 1's two nearest are {2,3} at
            % distances {h, 2h}, but vertex 3's two nearest are {2,4}, TIED at distance h
            % each -- vertex 1 (at distance 2h from vertex 3) does not make vertex 3's list.
            % So the 1-3 edge is one-way: it exists only from vertex 1's row. max(W,W') keeps
            % it at exp(-(2h)^2/(2h^2)) = exp(-2); averaging would halve it to exp(-2)/2.
            % This is what makes the test discriminate max from average -- plain symmetry
            % cannot, since both operations produce a symmetric matrix.
            G = rheome.sensors.graph(rheome.sensors.linear('NumSites', 8, 'Pitch', 1e-3), 'K', 2);
            W = full(G.W);
            tc.verifyEqual(W, W.', 'AbsTol', 1e-14);
            tc.verifyEqual(W(1,3), exp(-2), 'RelTol', 1e-12);
        end

        function defaultNeighbourCountFollowsIntrinsicDimension(tc)
            tc.verifyEqual(rheome.sensors.graph(rheome.sensors.linear('NumSites', 10)).K, 4);
            tc.verifyEqual(rheome.sensors.graph(rheome.sensors.grid('Size', [5 5], 'Pitch', 1e-3)).K, 8);
        end

        function defaultSigmaIsThePitch(tc)
            arr = rheome.sensors.grid('Size', [5 5], 'Pitch', 2e-3);
            tc.verifyEqual(rheome.sensors.graph(arr).Sigma, arr.Pitch, 'RelTol', 1e-12);
        end

        function surfaceShapeIsReturnedSoMeshCodeAcceptsIt(tc)
            % The whole point of the struct shape: existing readouts take it unchanged.
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [5 5], 'Pitch', 1e-3));
            tc.verifyTrue(all(isfield(G, {'Vertices','Faces','nV'})));
            tc.verifyEqual(G.nV, 25);
            tc.verifyEqual(G.Vertices, G.Array.Pos);
        end

        function normalizedLaplacianHasSpectrumInZeroTwo(tc)
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3), ...
                              'Laplacian', 'normalized');
            ev = eig(full(G.L));
            tc.verifyGreaterThan(min(ev), -1e-10);
            tc.verifyLessThan(max(ev), 2 + 1e-10);
        end

        function anUnknownLaplacianErrors(tc)
            tc.verifyError(@() rheome.sensors.graph(rheome.sensors.grid('utah'), 'Laplacian', 'nope'), ...
                'sensors:graph:laplacian');
        end

        function eigenvaluesAreAscendingAndStartAtZero(tc)
            M = rheome.sensors.modes(rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3)));
            tc.verifyEqual(M.Lambda, sort(M.Lambda), 'AbsTol', 1e-14);
            tc.verifyEqual(M.Lambda(1), 0, 'AbsTol', 1e-10);
        end

        function eigenvectorsAreOrthonormal(tc)
            M = rheome.sensors.modes(rheome.sensors.graph(rheome.sensors.grid('Size', [5 5], 'Pitch', 1e-3)));
            tc.verifyEqual(M.Phi.' * M.Phi, eye(25), 'AbsTol', 1e-10);
        end

        function theEigenTransformRoundTrips(tc)
            M = rheome.sensors.modes(rheome.sensors.graph(rheome.sensors.grid('Size', [5 5], 'Pitch', 1e-3)));
            X = randn(25, 3);
            tc.verifyEqual(M.T.inverse(M.T.forward(X)), X, 'AbsTol', 1e-9);
        end

        function chebyshevAndEigenAgreeOnAHeatFilter(tc)
            % The two routes must be the same operator. A heat kernel is the right probe:
            % smooth, so the polynomial converges, and it exercises the whole spectrum.
            G  = rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3));
            Me = rheome.sensors.modes(G, 'Route', 'eigen');
            Mc = rheome.sensors.modes(G, 'Route', 'chebyshev', 'Order', 80);
            X  = randn(36, 2);
            g  = @(l) exp(-0.4 * l / Me.Lmax);
            Ye = Me.Phi * (g(Me.Lambda) .* (Me.Phi.' * X));
            Yc = Mc.T.filter(g, X);
            tc.verifyEqual(Yc, Ye, 'RelTol', 1e-4, 'AbsTol', 1e-8);
        end

        function chebyshevLmaxBoundsTheTrueSpectrum(tc)
            % A Chebyshev lmax BELOW the operator's true lambda_max diverges rather than
            % losing accuracy, so this is a safety property, not an accuracy one. A >=
            % against the exact maximum is worthless here: both sides come from the same
            % deterministic eig on the same matrix, so exact equality would pass and the
            % margin -- the entire point of this test's name -- would go unguarded. Demand
            % a real margin instead, so a "tidy" of 1.01*max(Lam) down to max(Lam) fails.
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [6 6], 'Pitch', 1e-3));
            M = rheome.sensors.modes(G, 'Route', 'chebyshev');
            tc.verifyGreaterThan(M.Lmax, 1.005 * max(eig(full(G.L))));
        end

        function anUnknownRouteErrors(tc)
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [4 4], 'Pitch', 1e-3));
            tc.verifyError(@() rheome.sensors.modes(G, 'Route', 'nope'), 'sensors:modes:route');
        end

    end
end

% Author: Diellor Basha, 2026
