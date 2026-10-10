classdef tVortexSpectrum < matlab.unittest.TestCase
    % rheome.flow.vortexspectrum: matching a real sensor PSD with a background plus band vortices.
    %
    % Author: Diellor Basha, 2026

    properties
        B; g; st; dd; out
    end

    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                t.B  = rheome.load.bases(rheomeTestSubject());
                t.g  = rheome.operators.gauge(t.B.L.S.Vertices, double(t.B.L.S.Faces), Method="diffusion");
                t.st = rheome.load.study(rheomeTestSubject());
                t.dd = rheome.load.dirac(rheomeTestSubject());
            catch
                t.assumeFail('cached bases/study/dirac for test subject are not present');
            end
            t.out = rheome.flow.vortexspectrum(rheomeTestSubject(), Bands=[8 16; 16 32], Duration=8, ...
                        Bases=t.B, Gauge=t.g, Study=t.st, Dirac=t.dd, Verbose=false);
        end
    end

    methods (Test)

        function theBandTableIsWellFormed(t)
            b = t.out.bands;
            t.verifyEqual(height(b), 2);
            t.verifyTrue(all(b.used));
            t.verifyTrue(all(b.momentNAm > 0), 'every used band must carry a rhythm');
            t.verifyTrue(all(b.nAtoms >= 1));
            t.verifyEqual(b.fc(1), 8*sqrt(2), 'RelTol', 0.02);
            % constant-Q: an octave up halves the support
            t.verifyEqual(b.supportSec(1)/b.supportSec(2), 2, 'RelTol', 0.05);
        end

        function eachBandCarriesItsMeasuredPower(t)
            % ⭐ the amplitude is solved, not swept, so this should land close to one
            r = t.out.check.bandPowerRatio(t.out.bands.used);
            t.verifyEqual(r, ones(size(r)), 'AbsTol', 0.15);
        end

        function theSynthesisedExponentMatchesTheReal(t)
            t.verifyEqual(t.out.check.chiSim, t.out.check.chiReal, 'AbsTol', 0.4);
            t.verifyGreaterThan(t.out.check.logCorr, 0.6);
        end

        function theFitIsActuallyAFloor(t)
            % ⚠⚠ the whole synthesis depends on this: at rheome.spectral.aperiodic's DEFAULTS the fit
            %    sits above the data at 49.8% of bins and delta/theta get no rhythm at all
            t.verifyLessThan(t.out.check.fitOverFrac, 0.20, ...
                'the background must sit under the data, or there is no room for a rhythm');
            t.verifyGreaterThan(t.out.check.fitR2, 0.9);
        end

        function aBandTooSlowForTheRecordIsSkippedNotTruncated(t)
            % delta supports 5.05 s, so an 8 s record cannot carry it
            f = @() rheome.flow.vortexspectrum(rheomeTestSubject(), Bands=[2 4], Duration=8, ...
                        Bases=t.B, Gauge=t.g, Study=t.st, Dirac=t.dd, Verbose=false);
            t.verifyWarning(f, 'flow:vortexspectrum:shortRecord');
            w = warning('off', 'flow:vortexspectrum:shortRecord');
            c = onCleanup(@() warning(w));
            o = f();
            t.verifyFalse(o.bands.used(1));
            t.verifyTrue(isnan(o.bands.momentNAm(1)));
        end

        function theOutputRecordHasTheRightShape(t)
            t.verifySize(t.out.B, size(t.out.Breal));
            t.verifyTrue(all(isfinite(t.out.B(:))));
            t.verifyGreaterThan(norm(t.out.B, 'fro'), 0);
        end

    end
end

% Author: Diellor Basha, 2026
