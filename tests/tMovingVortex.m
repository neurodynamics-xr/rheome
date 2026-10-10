classdef tMovingVortex < matlab.unittest.TestCase
    % rheome.flow.movingvortex: a vortex core that travels, and what motion costs in band.
    %
    % Author: Diellor Basha, 2026

    properties
        B; g
    end

    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                t.B = rheome.load.bases(rheomeTestSubject());
                t.g = rheome.operators.gauge(t.B.L.S.Vertices, double(t.B.L.S.Faces), Method="diffusion");
            catch
                t.assumeFail('cached bases for test subject are not present');
            end
        end
    end

    methods (Test)

        function zeroSpeedDegeneratesToTheStaticAtom(t)
            m = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=0);
            a = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false);
            t.verifyEqual(m.check.nWaypoints, 1);
            t.verifyEqual(m.pathLengthMM, 0, 'AbsTol', 1e-12);
            t.verifyEqual(m.JA, a.JA, 'AbsTol', 1e-12);
            t.verifyEqual(m.check.inbandFrac, 1, 'AbsTol', 0.01);
            t.verifyEqual(m.check.passageSec, Inf, 'a still vortex never passes');
        end

        function theFieldIsExactlyRankTwoK(t)
            m = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=0.05);
            K = m.check.nWaypoints;
            t.verifyGreaterThan(K, 1);
            t.verifySize(m.JA, [size(m.JA,1) K]);
            t.verifySize(m.W,  [K m.nT]);
            J = m.materialise();
            for k = [1 round(m.nT/2) m.nT]
                t.verifyEqual(J(:,k), m.JA*(m.W(:,k)*m.a(k)) + m.JB*(m.W(:,k)*m.b(k)), ...
                              'AbsTol', 1e-15);
            end
        end

        function theWeightsArePartitionOfUnityWhileTravelling(t)
            m = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=0.05);
            cs = sum(m.W, 1);
            t.verifyEqual(cs, ones(1, m.nT), 'AbsTol', 1e-12, 'weights do not sum to one');
            t.verifyGreaterThanOrEqual(min(m.W(:)), 0, 'weights must not go negative');
        end

        function theCoreActuallyMoves(t)
            still = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=0);
            fast  = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=0.25);
            t.verifyEqual(still.check.speedEuclideanMS, 0, 'AbsTol', 1e-6);
            t.verifyGreaterThan(fast.check.speedEuclideanMS, 0.05);
            t.verifyGreaterThan(fast.pathLengthMM, 50);
        end

        function theGeodesicAndEuclideanSpeedsDifferByTheFolding(t)
            m = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=0.25);
            % ⚠ the core covers more SURFACE than it displaces through SPACE, so the Euclidean
            %   centroid speed must come in BELOW the geodesic one by the fold factor
            t.verifyGreaterThan(m.check.foldFactor, 1.2);
            t.verifyLessThan(m.check.speedEuclideanMS, m.check.speedGeodesicMS);
            t.verifyEqual(m.check.speedEuclideanMS*m.check.foldFactor, ...
                          m.check.speedGeodesicMS, 'RelTol', 0.35);
        end

        function speedBroadensTheBandMonotonically(t)
            % ⭐ the measurement the function exists to make
            v = [0.05 0.25 1.0];  frac = zeros(1,3);
            w = warning('off', 'flow:movingvortex:clamped');
            c = onCleanup(@() warning(w));
            for i = 1:3
                m = rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=v(i));
                frac(i) = m.check.inbandFrac;
            end
            t.verifyTrue(all(diff(frac) < 0), 'in-band energy must fall as speed rises');
            t.verifyGreaterThan(frac(1), 0.99, 'a slow vortex should stay in band');
            t.verifyLessThan(frac(3), 0.6, 'at 1 m/s most energy should have left the octave');
        end

        function tooFastToFitTheMeshWarnsAndClamps(t)
            f = @() rheome.flow.movingvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, SpeedMS=5);
            t.verifyWarning(f, 'flow:movingvortex:clamped');
            w = warning('off', 'flow:movingvortex:clamped');
            c = onCleanup(@() warning(w));
            m = f();
            t.verifyTrue(m.check.clamped);
            t.verifyLessThan(m.pathLengthMM, 200, 'the mesh is only ~161 mm across');
        end

    end
end

% Author: Diellor Basha, 2026
