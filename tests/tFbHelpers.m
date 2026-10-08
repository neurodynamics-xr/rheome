classdef tFbHelpers < matlab.unittest.TestCase
% The self-checking machinery the three demos rely on.

    methods (Test)

        function checkPassesWithinTolerance(tc)
            r = rheome.demos.fbd_selftest('check', 'width', 20.4, 20.0, 1.0, 'mm');
            tc.verifyTrue(r.pass);
            tc.verifyEqual(r.measured, 20.4);
            tc.verifyEqual(r.expected, 20.0);
        end

        function checkFailsOutsideTolerance(tc)
            r = rheome.demos.fbd_selftest('check', 'width', 25.0, 20.0, 1.0, 'mm');
            tc.verifyFalse(r.pass);
        end

        function checkPrintsAnAlignedLine(tc)
            s = evalc("rheome.demos.fbd_selftest('check','width',20.4,20.0,1.0,'mm');");
            tc.verifySubstring(s, 'width');
            tc.verifySubstring(s, 'PASS');
            tc.verifySubstring(s, 'mm');
        end

        function checkHandlesNaNAsFailure(tc)
            % A NaN measurement is a failure, never a silent pass.
            r = rheome.demos.fbd_selftest('check', 'slope', NaN, 10, 1, 'm/s');
            tc.verifyFalse(r.pass);
        end

        function reportSummarisesAndFlags(tc)
            good = rheome.demos.fbd_selftest('check', 'a', 1, 1, 0.1, '');
            bad  = rheome.demos.fbd_selftest('check', 'b', 5, 1, 0.1, '');
            [ok, n] = rheome.demos.fbd_selftest('report', [good good]);
            tc.verifyTrue(ok);   tc.verifyEqual(n, 2);
            [ok2, ~] = rheome.demos.fbd_selftest('report', [good bad]);
            tc.verifyFalse(ok2);
        end

        function ternaryPicksTheRightBranch(tc)
            tc.verifyEqual(rheome.demos.fbd_selftest('ternary', true,  'off', 'on'), 'off');
            tc.verifyEqual(rheome.demos.fbd_selftest('ternary', false, 'off', 'on'), 'on');
        end

        function figureIsInvisibleWhenExporting(tc)
            % Visible is a matlab.lang.OnOffSwitchState in R2023b, not a char.
            fH = rheome.demos.fbd_selftest('fig', true, [10 10 300 200]);
            c = onCleanup(@() close(fH));
            tc.verifyEqual(char(fH.Visible), 'off');
            fV = rheome.demos.fbd_selftest('fig', false, [10 10 300 200]);
            c2 = onCleanup(@() close(fV));
            tc.verifyEqual(char(fV.Visible), 'on');
        end

    end
end

% Author: Diellor Basha, 2026
