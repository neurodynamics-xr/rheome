classdef tForwardSimulate < matlab.unittest.TestCase
% rheome.forward.simulate: the mode-space forward model, and regressions for the two silent errors that
% produced complete, plausible, wrong tables before being caught.
%
% Author: Diellor Basha, 2026
    properties
        ctx
    end
    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject', 'noise');   % skips, with the reason, when the cache or the subject is absent
            try
                B = rheome.load.bases(rheomeTestSubject());
                t.ctx = struct('H', B.L, 'd', rheome.load.dirac(rheomeTestSubject()), ...
                    'st', rheome.load.study(rheomeTestSubject()), 'er', rheome.load.study(rheomeTestSubject('noise')));
            catch
                t.assumeFail('cached test subject / noise study data not present');
            end
        end
    end
    methods (Test)
        function modeSpaceAndVertexForwardsAgree(t)
            % ⭐ the internal check that catches using the wrong gain convention
            o = rheome.forward.simulate(rheomeTestSubject(), Duration=1, NoiseTrials=5, MomentNAm=1, ...
                Context=t.ctx, Verbose=false);
            t.verifyLessThan(o.routeRelErr, 1e-10, ...
                'B = Gmode*C must equal Gain*(Psi*C) exactly');
        end
        function theGainIsSynthesisNotAnalysis(t)
            % ⚠⚠ REGRESSION. rheome.forward.dirac is the ANALYSIS gain: it carries a mass weighting and
            % differs from the synthesis operator Gain*Phi by ~1/(vertex area) ~ 9e4. Substituting it
            % type-checks, runs, returns tesla and is 100+ dB low. This pins the distinction so the
            % two can never be quietly swapped.
            o = rheome.forward.simulate(rheomeTestSubject(), Duration=1, NoiseTrials=5, MomentNAm=1, ...
                Context=t.ctx, Verbose=false);
            % ⚠ the RAW gain includes bad and non-MEG channels whose rows are NaN, so
            % rheome.forward.dirac on it returns all-NaN and the comparison silently reads false.
            % rheome.forward.leadfield is the good-MEG restriction and exists for exactly this reason.
            Ga = rheome.forward.dirac(rheome.forward.leadfield(t.ctx.st), t.ctx.d);
            cols = find(t.ctx.d.Hemisphere == 1);
            % ⚠ omitnan: the near-null Dirac modes give NaN column norms in one of the two gains,
            % and a plain median then propagates it and the comparison reads NaN > 1000 = false.
            r = median(vecnorm(o.Gmode,2,1), 'omitnan') / median(vecnorm(Ga(:,cols),2,1), 'omitnan');
            t.verifyGreaterThan(r, 1e3, ...
                'the synthesis and analysis gains should differ by orders of magnitude');
        end
        function theSpatialScaleActuallyChangesTheField(t)
            % ⚠⚠ REGRESSION for the worst bug in this file's history: rheome.filters.mexhat was applied to
            % the DIRAC eigenvalues, which carry no length axis, so every requested wavelength gave
            % the same pattern and a 4x5 threshold table came back constant to three digits.
            a = rheome.forward.simulate(rheomeTestSubject(), WavelengthMM=267, Duration=1, NoiseTrials=5, ...
                MomentNAm=1, Context=t.ctx, Verbose=false);
            b = rheome.forward.simulate(rheomeTestSubject(), WavelengthMM=70, Duration=1, NoiseTrials=5, ...
                MomentNAm=1, Context=t.ctx, Verbose=false);
            t.verifyGreaterThan(abs(log10(a.rmsPerNAm/b.rmsPerNAm)), 0.05, ...
                'a 3.8x change in requested wavelength must change the sensor field');
        end
        function theThresholdIsClosedFormNotGridQuantised(t)
            % ⚠ reading the threshold off logspace(-1,2.5,15) quantises it to 1.78x steps, which
            % reported 3.2 nAm in all 20 cells of a band x scale sweep.
            o = rheome.forward.simulate(rheomeTestSubject(), Duration=1, NoiseTrials=20, ...
                MomentNAm=logspace(-1,2,13), Context=t.ctx, Verbose=false);
            t.verifyTrue(isfinite(o.threshMomentNAm) && o.threshMomentNAm > 0);
            t.verifyLessThan(abs(log2(o.threshMomentNAm/o.threshFromGrid)), 1.2, ...
                'the closed-form threshold must sit within one grid step of the swept one');
        end
        function theSensorFieldIsLinearInTheMoment(t)
            o = rheome.forward.simulate(rheomeTestSubject(), Duration=1, NoiseTrials=5, ...
                MomentNAm=[1 10], Context=t.ctx, Verbose=false);
            t.verifyEqual(o.rmsSignal(2)/o.rmsSignal(1), 10, 'RelTol', 1e-9);
        end
        function everyFamilyProducesAField(t)
            for f = ["oscillator" "resonator" "gabor" "travwave" "wave" "diffusion"]
                o = rheome.forward.simulate(rheomeTestSubject(), Family=f, Duration=1, NoiseTrials=5, ...
                    MomentNAm=1, Context=t.ctx, Verbose=false);
                t.verifyGreaterThan(o.rmsPerNAm, 0, f);
                t.verifyEqual(size(o.B,2), o.nT, f);
            end
        end
    end
end
