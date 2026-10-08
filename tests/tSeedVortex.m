classdef tSeedVortex < matlab.unittest.TestCase
% rheome.flow.seedvortex: imposing circulation in the connection frame. The two dials are the point.
%
% Author: Diellor Basha, 2026
    properties
        B; g
    end
    methods (TestClassSetup)
        function build(t)
            try
                t.B = rheome.load.bases(rheomeTestSubject());
                t.g = rheome.operators.gauge(t.B.L.S.Vertices, double(t.B.L.S.Faces), Method="diffusion");
            catch
                t.assumeFail('cached bases for test subject are not present');
            end
        end
    end
    methods (Test)
        function phaseSetsTheTypeAtFixedWinding(t)
            % ⭐⭐ THE WHOLE POINT, and what the stream-function route could not do. Winding +1 in both
            % cases; only the 90-degree phase offset decides source against vortex.
            a = rheome.flow.seedvortex(rheomeTestSubject(), Phase=0,    Bases=t.B, Gauge=t.g);
            b = rheome.flow.seedvortex(rheomeTestSubject(), Phase=pi/2, Bases=t.B, Gauge=t.g);
            t.verifyEqual(a.check.charge, 1);
            t.verifyEqual(b.check.charge, 1);
            t.verifyEqual(a.check.type, "source");
            t.verifyEqual(b.check.type, "vortex");
            t.verifyGreaterThan(b.check.curlOverDiv, 2*a.check.curlOverDiv);
        end
        function windingSetsTheSign(t)
            a = rheome.flow.seedvortex(rheomeTestSubject(), Winding=-1, Bases=t.B, Gauge=t.g);
            b = rheome.flow.seedvortex(rheomeTestSubject(), Winding=+1, Bases=t.B, Gauge=t.g);
            t.verifyEqual(a.check.charge, -1);
            t.verifyEqual(b.check.charge, +1);
            t.verifyEqual(a.check.type, "saddle");
        end
        function poincareHopfHolds(t)
            o = rheome.flow.seedvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g);
            t.verifyEqual(o.check.chi, 2);
        end
        function theFieldIsPurelyTangential(t)
            o = rheome.flow.seedvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g);
            t.verifyEqual(o.check.tangentialFrac, 1, 'AbsTol', 1e-9);
        end
        function theDefaultSeedAvoidsTheGaugeSingularities(t)
            % ⚠ round(nV/2) landed 8 mm from one, where the observed charge is m plus the gauge's own
            o = rheome.flow.seedvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g);
            t.verifyGreaterThan(o.check.gaugeSingularDistMM, 30);
        end
        function aFinerVortexIsCleanerThanACoarseOne(t)
            % ⚠ the opposite of what the instrument prefers: 2*r50 = 104 mm, but the coarse pattern
            % spans more curvature and fragments into more strong critical points
            a = rheome.flow.seedvortex(rheomeTestSubject(), WavelengthMM=70,  Bases=t.B, Gauge=t.g);
            b = rheome.flow.seedvortex(rheomeTestSubject(), WavelengthMM=267, Bases=t.B, Gauge=t.g);
            t.verifyLessThan(a.check.nStrong, b.check.nStrong);
            t.verifyLessThanOrEqual(a.check.nStrong, 2);
        end
        function curlOverDivStaysNearOneAndThatIsTheSurface(t)
            % ⚠⚠ REGRESSION FOR A CLAIM I WITHDREW. Imposing the azimuthal direction does NOT avoid the
            % curvature divergence: a unit azimuthal field on a curved surface has divergence from the
            % curvature itself, so ~1.1-1.2 is near the floor and no construction beats it here.
            o = rheome.flow.seedvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g);
            t.verifyGreaterThan(o.check.curlOverDiv, 0.9);
            t.verifyLessThan(o.check.curlOverDiv, 3, ...
                'if this ever exceeds 3 the operators changed and the docstring is stale');
        end
    end
end
