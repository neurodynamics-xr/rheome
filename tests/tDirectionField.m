classdef tDirectionField < matlab.unittest.TestCase
% rheome.flow.directionfield: planted singularities, verified by holonomy AND by the field's own winding.
%
% Author: Diellor Basha, 2026
    properties
        B; C
    end
    methods (TestClassSetup)
        function build(t)
            try
                t.B = rheome.load.bases(rheomeTestSubject());
                t.C = rheome.operators.connection_laplacian(t.B.L.S.Vertices, double(t.B.L.S.Faces));
            catch
                t.assumeFail('cached bases for test subject are not present');
            end
        end
    end
    methods (Test)
        function theHolonomyIsExactlyWhatWasPrescribed(t)
            % ⭐⭐ THE CHECK THAT MATTERS for a designed field: not detection, holonomy.
            o = rheome.flow.directionfield(rheomeTestSubject(), Bases=t.B, Connection=t.C);
            t.verifyEqual(o.check.holonomyAtPrescribed, double(o.charges), 'AbsTol', 1e-6);
            t.verifyLessThan(o.check.holonomyElsewhereMax, 1e-9, ...
                'the connection must be flat away from the prescribed singularities');
        end
        function theFieldIsExactlyParallelAwayFromThem(t)
            o = rheome.flow.directionfield(rheomeTestSubject(), Bases=t.B, Connection=t.C);
            t.verifyLessThan(o.check.transportResidual, 1e-8);
        end
        function chargesMustSumToChi(t)
            % ⚠ the wrapped face curvatures sum to 2*pi*chi; [1 -1] was once required and is now wrong
            t.verifyError(@() rheome.flow.directionfield(rheomeTestSubject(), Faces=[10 20], Charges=[1 -1], ...
                Bases=t.B, Connection=t.C), 'flow:directionfield:charges');
        end
        function theFieldIsPurelyTangentialAndUnitBeforeEnveloping(t)
            o = rheome.flow.directionfield(rheomeTestSubject(), Bases=t.B, Connection=t.C);
            t.verifyEqual(o.check.tangentialFrac, 1, 'AbsTol', 1e-9);
            t.verifyEqual(abs(o.z), ones(size(o.z)), 'AbsTol', 1e-9);
        end
        function rotatingDoesNotMoveTheSingularities(t)
            % ⭐ Rotate adds a constant to phi, so the holonomy is untouched: the 90-degree turn
            % changes what the field IS, never where its defects are.
            a = rheome.flow.directionfield(rheomeTestSubject(), Rotate=0,    Bases=t.B, Connection=t.C);
            b = rheome.flow.directionfield(rheomeTestSubject(), Rotate=pi/2, Bases=t.B, Connection=t.C);
            t.verifyEqual(a.check.holonomyAtPrescribed, b.check.holonomyAtPrescribed, 'AbsTol', 1e-9);
            t.verifyEqual(a.faces, b.faces);
        end
        function theFieldWindsExactlyWherePrescribed(t)
            % ⭐⭐ REGRESSION FOR THE UNWRAPPED SOLVE. This test used to pin the opposite -- "the per-face
            % winding is spread over >100 faces, a detector limitation" -- and that was the bug: the
            % index was k - n_f. With the wrap the field's own winding is the prescription, exactly.
            o = rheome.flow.directionfield(rheomeTestSubject(), Bases=t.B, Connection=t.C);
            t.verifyEqual(o.check.windingSumAll, 2);
            t.verifyEqual(o.check.windingNonzeroCount, numel(o.faces));
            t.verifyEqual(o.check.windingAtPrescribed, double(o.charges));
        end
        function aVortexAntivortexPairBesideTheBudget(t)
            o = rheome.flow.directionfield(rheomeTestSubject(), Faces=[10 5000 12000 18000], Charges=[1 -1 1 1], ...
                Bases=t.B, Connection=t.C);
            t.verifyEqual(o.check.windingAtPrescribed, [1 -1 1 1]);
            t.verifyEqual(o.check.windingNonzeroCount, 4);
        end
    end
end
