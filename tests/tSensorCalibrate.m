classdef tSensorCalibrate < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (Test)

        function sampleEvaluatesAFieldAtEverySensorAndFrame(tc)
            arr = rheome.sensors.grid('Size', [3 3], 'Pitch', 1e-3);
            U   = rheome.sensors.sample(arr, @(P, tj) P(:,1) + tj, [0 1 2]);
            tc.verifyEqual(size(U), [9 3]);
            tc.verifyEqual(U(:,2), arr.Pos(:,1) + 1, 'AbsTol', 1e-15);
        end

        function sampleWithoutATimeAxisReturnsOneFrame(tc)
            arr = rheome.sensors.grid('Size', [3 3], 'Pitch', 1e-3);
            U   = rheome.sensors.sample(arr, @(P) P(:,2));
            tc.verifyEqual(size(U), [9 1]);
            tc.verifyEqual(U, arr.Pos(:,2), 'AbsTol', 1e-15);
        end

        function closedFormAlphaMatchesTheHandComputedInteriorValue(tc)
            % An interior sensor of a square lattice with K=8 and sigma=h keeps 4 edge
            % neighbours at h (w = e^-1/2) and 4 diagonals at h*sqrt(2) (w = e^-1). So
            %   sum_j w_ij |d_ij|^2 = h^2 (4 e^-1/2 + 8 e^-1)
            % and alpha = that / (2*dim) with dim = 2. Hand-computable, so this is a real
            % test of the implementation rather than a self-consistency check.
            h   = 1e-3;
            arr = rheome.sensors.grid('Size', [11 11], 'Pitch', h);
            G   = rheome.sensors.graph(arr, 'K', 8, 'Sigma', h);
            c   = rheome.sensors.calibrate(G);
            want = (h^2 * (4*exp(-0.5) + 8*exp(-1))) / 4;
            interior = sub2ind([11 11], 6, 6);
            tc.verifyEqual(c.alphaPerVertex(interior), want, 'RelTol', 1e-12);
        end

        function theTwoRoutesAgreeOnALargeRegularLattice(tc)
            % The closed form assumes every sensor sees the same neighbourhood, which is
            % true only in the interior. The two routes therefore converge as the boundary
            % fraction shrinks -- hence 32x32 (12% boundary), not 6x6 (89%).
            h = 1e-3;
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [32 32], 'Pitch', h), 'K', 8, 'Sigma', h);
            c = rheome.sensors.calibrate(G);
            tc.verifyLessThan(c.agree, 0.05);
        end

        function agreementDegradesOnASmallArrayAndThatIsTheBoundary(tc)
            % Spec risk 2, pinned rather than hoped about. A 12x12 lattice is far more
            % boundary than a 32x32 one, so its two routes agree less well -- and a test
            % that asserted otherwise would be asserting the boundary away.
            %
            % 12x12, not 10x10: a 10x10 lattice's Aperture/Pitch (12.7) falls just short of
            % a small-k fit window and silently hits calibrate's 3-point fallback (see
            % sensors:calibrate:narrowFitWindow), which has its OWN bias on top of the
            % boundary effect this test is about. 12x12 clears the window (nFit=4) with
            % room to spare, so both .hasFitWindow flags below being true is what keeps
            % this comparison about the boundary and nothing else.
            h  = 1e-3;
            cB = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [32 32], 'Pitch', h)));
            cS = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [12 12], 'Pitch', h)));
            tc.verifyTrue(cB.hasFitWindow);
            tc.verifyTrue(cS.hasFitWindow);
            tc.verifyGreaterThan(cS.agree, cB.agree);
            tc.verifyLessThan(cS.agree, 0.15);
        end

        function lambdaUsableIsCoarserThanTheAliasFloor(tc)
            % THE ALIAS FLOOR IS NOT THE ACCURACY FLOOR. lambda ~ alpha k^2 is a small-k
            % expansion; it is already 10% low around 5-6 sensors per cycle, well before
            % the 2-sensor alias limit. lambdaUsable is where the model breaks, and it is
            % the number a wavelength claim actually rests on.
            h = 1e-3;
            c = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [32 32], 'Pitch', h), ...
                                                'K', 8, 'Sigma', h));
            tc.verifyGreaterThan(c.lambdaUsable, 2*h);      % coarser than the alias floor
            tc.verifyGreaterThan(c.lambdaUsable, 4*h);      % several sensors per cycle
            tc.verifyLessThan(c.lambdaUsable, 12*h);
        end

        function lambdaUsableOnAChainIsSustainedNotAnIsolatedBlip(tc)
            % A 1D chain has exactly one measurement direction (see i_directions), so
            % lambda_hat(k) has none of the 8-direction averaging that damps route 2's
            % finite-domain artifact on a planar array. On a 64-site probe that curve dips
            % back under 10% repeatedly near kMin (isolated blips at ~10.05% and ~10.24%)
            % before turning over for good at ~12.05% and climbing monotonically -- the
            % FIRST crossing alone would land lambdaUsable at ~37 sensors/cycle, five times
            % coarser than this file's own header (":26-30") documents as "10% low at
            % roughly 5-6 sensors per cycle". Requiring the crossing to hold at the next
            % grid point too is what this test guards: nothing else in this suite exercises
            % a Dim-1 array's departure curve.
            arr = rheome.sensors.linear('NumSites', 64, 'Pitch', 100e-6);
            c   = rheome.sensors.calibrate(rheome.sensors.graph(arr));
            spc = c.lambdaUsable / arr.Pitch;         % sensors per cycle
            tc.verifyGreaterThan(spc, 5);
            tc.verifyLessThan(spc, 12);
        end

        function alphaGivesBackThePlantedWavelength(tc)
            % The claim the whole phase rests on: a planted wave of known wavelength must
            % read back as that wavelength, in metres, through lambda_hat and alpha alone.
            %
            % Tolerance is 0.08, not 0.06. The lattice model gives
            % 2(1-cos(kh))/(kh)^2 = 0.9497 at 8 sensors per cycle, so the recovered
            % wavelength is already ~2.6% high from the quadratic model's OWN curvature --
            % not noise -- before the boundary contamination that a finite 24x24 array puts
            % into the fitted alpha is added on top. 0.06 leaves almost no headroom over
            % that known systematic.
            h  = 1e-3;
            G  = rheome.sensors.graph(rheome.sensors.grid('Size', [24 24], 'Pitch', h), 'K', 8, 'Sigma', h);
            c  = rheome.sensors.calibrate(G);
            lamWant = 8e-3;                                  % 8 mm, 8 sensors per cycle
            k  = 2*pi/lamWant;
            u  = rheome.sensors.sample(G.Array, @(P) cos(k * P(:,1)));
            u  = u - mean(u);
            lh = (u.' * (G.L * u)) / (u.' * u);
            lamGot = 2*pi / sqrt(lh / c.alpha);
            tc.verifyEqual(lamGot, lamWant, 'RelTol', 0.08);
        end

        function alphaSpreadTracksBoundaryFractionNotJustIrregularity(tc)
            % alphaSpread is neighbourhood heterogeneity, and at these array sizes BOUNDARY
            % FRACTION is the dominant term -- not "regular lattice vs irregular cap" as a
            % story on its own. A 32x32 lattice (12% boundary) has markedly lower spread
            % than a 10x10 lattice (36% boundary); and the 64-point EEG cap, despite being
            % geometrically irregular, has a SMALLER boundary fraction than the 10x10 grid
            % and a correspondingly smaller alphaSpread than it -- the opposite of what a
            % naive "lattice is regular, cap is irregular" framing would predict. (An
            % earlier version of this test compared the cap only to a 16x16 lattice, which
            % happens to have a still-larger spread than the cap; that comparison doesn't
            % generalise, so it's replaced with the boundary-fraction story the numbers
            % actually support.
            %
            % Both the 10x10 lattice and the 64-channel cap are below route 2's small-k
            % fit window (see sensors:calibrate:narrowFitWindow) and warn accordingly; that
            % warning is about .alpha, not .alphaSpread (a route-1-only quantity), so it is
            % suppressed here rather than being this test's point.
            h  = 1e-3;
            w  = warning('off', 'sensors:calibrate:narrowFitWindow');
            oc = onCleanup(@() warning(w));
            cBig   = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [32 32], 'Pitch', h)));
            cSmall = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [10 10], 'Pitch', h)));
            cCap   = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.eeg('NumChannels', 64)));
            tc.verifyLessThan(cBig.alphaSpread, cSmall.alphaSpread);
            tc.verifyLessThan(cCap.alphaSpread, cSmall.alphaSpread);
        end

        function speedScaleIsTheSquareRootOfAlpha(tc)
            % This pins the assignment (speedScale = sqrt(alpha)) and catches a typo on
            % that line -- it does NOT pin the direction of the physical-speed rule
            % (dispersion() fits omega = c_graph*sqrt(lambda); since sqrt(lambda) =
            % sqrt(alpha)*k, physical speed is c_graph*speedScale, a MULTIPLICATION, per
            % the header). An end-to-end check -- plant a wave of known speed, run
            % rheome.dynamics.dispersion, and confirm c_graph*speedScale recovers it -- belongs to
            % a later task once rheome.dynamics.dispersion is calibrated against this operator.
            c = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [12 12], 'Pitch', 1e-3)));
            tc.verifyEqual(c.speedScale, sqrt(c.alpha), 'RelTol', 1e-12);
        end

        function calibrateRefusesTheNormalizedLaplacian(tc)
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [8 8], 'Pitch', 1e-3), ...
                              'Laplacian', 'normalized');
            tc.verifyError(@() rheome.sensors.calibrate(G), 'sensors:calibrate:laplacian');
        end

        function theDispersionCurveDepartsFromAlphaKSquaredNearTheAliasFloor(tc)
            % lambda_hat(k) tracks alpha*k^2 at long wavelengths and falls BELOW it as the
            % wavelength approaches 2*pitch. That departure is the resolution floor, and it
            % is derived here rather than asserted.
            h = 1e-3;
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [24 24], 'Pitch', h), 'K', 8, 'Sigma', h);
            c = rheome.sensors.calibrate(G);
            model = c.alpha * c.kGrid.^2;
            ratio = c.lambdaHat ./ model;
            tc.verifyEqual(ratio(1), 1, 'AbsTol', 0.1);          % long wavelength: on model
            tc.verifyLessThan(ratio(end), 0.85);                 % near the floor: below it
        end

        function hasFitWindowIsFalseOnACoarsePresetAndTrueOnALargeLattice(tc)
            % ECoG (8x8 @ 10 mm) is 9.9 pitches across its own diagonal: no wavenumber can
            % be both >=1.5 cycles across that aperture (route 2's own long-wavelength
            % margin) and inside the small-k fit window, so route 2 has NO small-k regime
            % on this array at all -- and calibrate must say so (warn + .hasFitWindow)
            % rather than answer as though a fit window existed.
            Gecog = rheome.sensors.graph(rheome.sensors.grid('ecog'));
            tc.verifyWarning(@() rheome.sensors.calibrate(Gecog), 'sensors:calibrate:narrowFitWindow');

            w = warning('off', 'sensors:calibrate:narrowFitWindow');
            oc = onCleanup(@() warning(w));
            cEcog = rheome.sensors.calibrate(Gecog);
            tc.verifyFalse(cEcog.hasFitWindow);
            % The fallback's actual fit window sits well past the nominal request --
            % .kFit (what was really fitted) and .kFitRequested must disagree here.
            tc.verifyGreaterThan(cEcog.kFit, cEcog.kFitRequested);

            h    = 1e-3;
            cBig = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('Size', [32 32], 'Pitch', h)));
            tc.verifyTrue(cBig.hasFitWindow);
            tc.verifyGreaterThanOrEqual(cBig.nFit, 3);
        end

        function calibrateRefusesAnApertureBelowThreeTimesThePitch(tc)
            % kMin (1.5 cycles across the aperture) exceeds kMax (the pi/Pitch alias floor)
            % whenever Aperture < 3*Pitch -- logspace would run backwards and the fallback
            % would fit the three WORST points on the grid. A 2x2 lattice (Aperture/Pitch =
            % sqrt(2) ~ 1.41) is exactly this case, and rheome.sensors.grid accepts it.
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [2 2], 'Pitch', 1e-3));
            tc.verifyError(@() rheome.sensors.calibrate(G), 'sensors:calibrate:aperture');
        end

        function lambdaUsableIsNaNWhenNoUsableBandExists(tc)
            % NaN requires the departure to be SUSTAINED across the two longest-wavelength
            % probed points, not merely crossed once at the first -- a coarse wavenumber
            % grid (NumK below the default) on the already-fallback ECoG preset makes both
            % point 1 (ratio 1.118) AND point 2 (ratio 1.103) exceed 10%, which is what
            % actually triggers NaN here, not point 1 alone.
            w  = warning('off', 'sensors:calibrate:narrowFitWindow');
            oc = onCleanup(@() warning(w));
            c  = rheome.sensors.calibrate(rheome.sensors.graph(rheome.sensors.grid('ecog')), 'NumK', 6);
            tc.verifyTrue(isnan(c.kUsable));
            tc.verifyTrue(isnan(c.lambdaUsable));
        end

    end
end

% Author: Diellor Basha, 2026
