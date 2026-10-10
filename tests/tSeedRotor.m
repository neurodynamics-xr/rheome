classdef tSeedRotor < matlab.unittest.TestCase
% rheome.flow.seedrotor: what a stream-function vortex actually guarantees. The topology is exact; the
% differential identities are discretisation-limited, and the tests say which is which.
%
% Author: Diellor Basha, 2026
    properties
        B
    end
    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try, t.B = rheome.load.bases(rheomeTestSubject());
            catch, t.assumeFail('cached bases for test subject are not present'); end
        end
    end
    methods (Test)
        function itIsNotAnIsolatedVortex(t)
            % ⚠⚠ .check reports the NEAREST critical point, which flattered the construction. At
            % 140 mm there are ~60 of them and the strongest is a saddle 8 mm from the seed, stronger
            % than the central vortex. Pinned so "a clean single rotor" is never claimed again.
            o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B);
            t.verifyGreaterThan(o.check.nCritical, 10, ...
                'if this ever drops to ~1 the operators changed and the docstring is stale');
        end
        function theWindingIsExactlyOneAtTheSeed(t)
            % ⭐ THE GUARANTEE THAT HOLDS. A maximum of the stream function is a centre, so the
            % rotated gradient has winding +1 -- and rheome.detect.criticalPoints gets winding exactly
            % (Poincare-Hopf), which is why this is sound ground truth for DETECTION.
            for wl = [267 140 70]
                o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=wl, Bases=t.B);
                t.verifyEqual(o.check.charge, 1, sprintf('%g mm', wl));
                t.verifyEqual(o.check.type, "vortex", sprintf('%g mm', wl));
                t.verifyLessThan(o.check.distMM, 10, 'the vortex must land on the seed');
            end
        end
        function theFieldIsNotAGradient(t)
            % the point of the cross product: grad(psi) would be 100% irrotational
            o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B);
            t.verifyLessThan(o.check.irrotationalFrac, 0.02);
            t.verifyGreaterThan(o.check.solenoidalFrac, 0.6);
        end
        function theTangentialProjectionRemovesTheNormalPartExactly(t)
            o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B, Tangential=true);
            t.verifyLessThan(o.check.normalFrac, 1e-9);
        end
        function theDivergenceIsNotZeroAndThatIsDiscretisation(t)
            % ⚠⚠ REGRESSION FOR A CLAIM I WITHDREW. A co-gradient is divergence-free intrinsically
            % and on a flat sheet; the ambient strong divergence of this field is ~1.3x the curl,
            % because the face-to-vertex average crosses folds and divergence carries curvature.
            % Pinned so nobody re-asserts div = 0 from the construction alone.
            o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B);
            t.verifyGreaterThan(o.check.divRatio, 0.5, ...
                'if this ever drops far below 1, the operators changed and the docstring is stale');
        end
        function theHelmholtzBudgetDoesNotClose(t)
            % ⚠ the parts are orthogonal only to ~20% discretely, so energies do not add
            o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=140, Bases=t.B);
            c = o.check;
            tot = c.solenoidalFrac + c.irrotationalFrac + c.normalFrac + c.harmFrac;
            t.verifyLessThan(tot, 0.95, 'the budget is not expected to close on this mesh');
            t.verifyGreaterThan(tot, 0.6);
        end
        function curlTracksMinusLaplacianBetterAtFineScales(t)
            % ⚠ corr 0.84 at 70 mm against 0.43 at 267 mm: a coarse pattern spans more curvature
            a = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=70,  Bases=t.B);
            b = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=267, Bases=t.B);
            t.verifyGreaterThan(a.check.curlCorr, b.check.curlCorr);
            t.verifyGreaterThan(a.check.curlCorr, 0.7);
        end
        function aLowPassKernelIsAllowedButNotLocalised(t)
            % ⚠ heat is low-pass, so its "vortex" is the whole hemisphere; allowed for comparison
            o = rheome.flow.seedrotor(rheomeTestSubject(), WavelengthMM=267, Kernel="heat", Bases=t.B);
            t.verifyEqual(o.check.charge, 1);
            t.verifyLessThan(o.check.curlCorr, 0.6);
        end
    end
end
