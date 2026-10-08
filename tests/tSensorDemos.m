classdef tSensorDemos < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function graphDemoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_graph();
            tc.verifyTrue(out.ok, 'every self-check in the sensor graph demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function graphDemoOpensFiveFigures(tc)
            close all;
            rheome.demos.sensor_graph();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function graphDemoAgreesOnAlphaForEveryRegularArray(tc)
            % ⚠ ONLY ARRAYS WITH A SMALL-K FIT WINDOW. The regression route needs >=3
            % wavenumbers below FitFraction*kMax, which exists only when the array spans
            % enough pitches. Utah and ECoG do NOT -- an 8x8 grid is 9.9 pitches across its
            % diagonal, so no wavenumber is both long enough to clear the aperture artefact
            % and long enough to clear the sampling floor. On those arrays .agree compares a
            % good estimate against a known-biased one and asserting it would be theatre.
            % The regular lattices are where the closed form is exact IN THE INTERIOR, so
            % they are where a disagreement is the code's fault rather than the geometry's.
            % 0.20 and not 0.05 because these arrays are small: tSensorCalibrate holds a
            % 32x32 lattice to 0.05, and the gap between the two numbers is the boundary
            % plus the fit window's sensitivity on a coarse array.
            close all;
            out = rheome.demos.sensor_graph();
            ok = out.regular & out.hasWindow;
            tc.verifyTrue(any(ok), 'at least one demo array must HAVE a small-k fit window');
            tc.verifyLessThan(max(out.alphaAgree(ok)), 0.20);
            % An array without the window is not a failure to hide -- assert it is REPORTED.
            tc.verifyFalse(all(out.hasWindow), ...
                'the demo must include at least one array with NO fit window, as a result');
        end

        function graphDemoShowsChebyshevAndEigenAgreeing(tc)
            close all;
            out = rheome.demos.sensor_graph();
            tc.verifyLessThan(out.chebErr, 1e-3);
        end

        function graphDemoDropsNoArrayFromTheSelfCheckAccounting(tc)
            % out.handled is populated BY the demo's three self-check loops, at the point
            % each array is actually dealt with -- NOT recomputed here from
            % regular/hasWindow. A test that recomputed the loop conditions would be a
            % tautology over two booleans (every element of a two-flag truth table falls
            % into exactly one of the three recomputed terms BY CONSTRUCTION, regardless of
            % what the loops do) and would still pass even if a loop were deleted. This is
            % exactly what happened to the MEG helmet before the loop existed: it fell
            % between the scored and no-window loops and nothing caught it.
            close all;
            out = rheome.demos.sensor_graph();
            tc.verifyTrue(all(out.handled), ...
                'every array must be visited by a scoring or reporting loop -- none silently dropped');
            tc.verifyNumElements(out.handled, numel(out.arrays));
        end

        function theDemoPrintsTheDisagreementForAnIrregularWindowedArray(tc)
            % Nothing else asserts the reporting fprintf's actually produce text, so a
            % refactor could turn them into no-ops silently. Capture stdout and check the
            % substance for the one array class (irregular, but WITH a fit window) whose
            % reporting loop was the whole point of this fix round.
            close all;
            txt = evalc('out = rheome.demos.sensor_graph();');
            irreg = find(~out.regular & out.hasWindow);
            if isempty(irreg)
                % No such array in this checkout (the MEG cache is absent) -- nothing to
                % assert. Returning early is honest here; asserting anyway would be a
                % vacuous pass dressed up as a real one.
                return;
            end
            tc.verifyNotEmpty(regexp(txt, 'DISAGREE', 'once'), ...
                'an irregular array with a good fit window must have its disagreement PRINTED');
            tc.verifyNotEmpty(regexp(txt, regexptranslate('escape', out.arrays{irreg(1)}.Name), 'once'));
        end

        function graphDemoMakesNoSourceClaim(tc)
            % Spec section 1.1 is load-bearing. A source word in this file is a defect.
            src = fileread(which('rheome.demos.sensor_graph'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|inverse|dipole|anatom', 'once'));
        end

        function limitsDemoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_limits();
            tc.verifyTrue(out.ok, 'every self-check in the sensor limits demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function limitsDemoOpensFourFigures(tc)
            close all;
            rheome.demos.sensor_limits();
            tc.verifyEqual(numel(findobj('Type','figure')), 4);
        end

        function insideTheWindowThePlantedWavelengthComesBack(tc)
            close all;
            out = rheome.demos.sensor_limits();
            tc.verifyTrue(any(arrayfun(@(s) any(s.inWindow), out.sweep)), ...
                'at least one array must HAVE a usable band, or this test is vacuous');
            for s = out.sweep
                sel = s.inWindow;
                if ~any(sel), continue; end   % no usable band -- reported, not swept
                tc.verifyLessThan(max(abs(s.lamRecovered(sel)./s.lamPlanted(sel) - 1)), 0.15, ...
                    sprintf('%s: inside its own window the array must recover the wavelength', s.name));
            end
        end

        function belowTheUsableFloorRecoveryFails(tc)
            % The floor must BITE. If short wavelengths came back correctly the floor would
            % be a decoration rather than a bound.
            close all;
            out = rheome.demos.sensor_limits();
            for s = out.sweep
                if ~any(s.inWindow), continue; end
                % 0.35, not 0.5: at 0.5 the mildest below-floor point carries only ~24%
                % error against a 20% threshold -- too close to call.
                below = s.lamPlanted < 0.35 * min(s.lamPlanted(s.inWindow));
                if ~any(below), continue; end
                r = s.lamRecovered(below);  q = s.lamPlanted(below);
                % Failure is EITHER a wrong number OR no number at all: at the alias floor
                % the field goes constant across the array and there is nothing to read.
                bad = ~isfinite(r) | abs(r./q - 1) > 0.2;
                tc.verifyTrue(all(bad), ...
                    sprintf('%s: below its usable floor recovery must fail', s.name));
            end
        end

        function theProbeReadsAnObliqueWaveAsCosThetaTimesK(tc)
            % A 1D array sees only k along its axis, so a wave crossing at theta reads
            % k*cos(theta) -- a LONGER apparent wavelength and a FASTER apparent speed.
            close all;
            out = rheome.demos.sensor_limits();
            tc.verifyEqual(out.obliqueK ./ out.obliqueK(1), cos(out.obliqueTheta), ...
                'RelTol', 0.05);
        end

        function limitsDemoMakesNoSourceClaim(tc)
            src = fileread(which('rheome.demos.sensor_limits'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|inverse|dipole|anatom', 'once'));
        end

    end
end

% Author: Diellor Basha, 2026
