classdef tJfbDispersion < matlab.unittest.TestCase

    methods (Test)

        function recoversTheSpeedOfASyntheticDiagonal(tc)
            % A structure travelling at c satisfies omega = c*sqrt(lambda), so it is a
            % DIAGONAL RIDGE in (sqrt(lambda), omega) and its slope IS the speed.
            % c must be chosen so the diagonal actually SPANS the frequency axis. On this
            % fixture k reaches 67 rad/m and omega reaches 804 rad/s, so c = 10 puts the
            % ridge at omega = 670 at the top wavenumber -- most of the axis. At c = 0.8
            % the ridge would reach only omega = 54, the bottom 7%, and there is no
            % diagonal to fit. That is the .fSupport caveat, not a defect in the fit.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',256, 'SamplingFrequency',256);
            lam = j.Lambda;  om = j.Omega;
            c = 10;
            C = exp(-(om - c*sqrt(lam)).^2 / (2*(0.03*max(om))^2));
            d = dispersion(j, C);
            tc.verifyEqual(d.slope, c, 'RelTol', 0.05);
            tc.verifyGreaterThan(d.slopeR2, 0.9);
        end

        function reportsTheSupportSoAClippedFitIsVisible(tc)
            % ⚠ THE BAND MUST SPAN THE DIAGONAL. A narrow retention sees a slice of a
            % diagonal spanning tens of Hz, and the slope is then not recoverable --
            % .fSupport is what makes that checkable instead of silent.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'Frequencies', linspace(8, 13, 24));
            lam = j.Lambda;  om = j.Omega;
            C = exp(-(om - 10*sqrt(lam)).^2 / (2*(0.03*max(om))^2));
            d = dispersion(j, C);
            tc.verifySize(d.fSupport, [1 2]);
            tc.verifyLessThanOrEqual(d.fSupport(2) - d.fSupport(1), 5);
        end

        function weightingIgnoresEmptyBins(tc)
            % Unweighted, bins carrying no signal dominate by sheer number and bias the
            % slope badly. Adding empty high-frequency bins must not move the answer much.
            fx = jfbFixture();
            jN = rheome.jointfilterbank(fx.gfb, 'SignalLength',256, 'SamplingFrequency',256);
            lam = jN.Lambda;  om = jN.Omega;
            C = exp(-(om - 10*sqrt(lam)).^2 / (2*(0.03*max(om))^2));
            C(:, jN.Frequencies > 115) = 0;                % explicitly empty bins
            d = dispersion(jN, C);
            tc.verifyEqual(d.slope, 10, 'RelTol', 0.1);
        end

        function degenerateInputReturnsNaNsNotAnError(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            d  = dispersion(j, zeros(numel(j.Lambda), j.NumOmega));
            tc.verifyTrue(isnan(d.slope));
            tc.verifyEmpty(d.kbar);
        end

        function anOffAxisRidgeReturnsNaNNotAZeroSlope(tc)
            % With the ridge almost entirely off the retained axis every frequency
            % reports the same mean wavenumber, so there is nothing to fit. A rank-
            % deficient least-squares would hand back slope = 0, which looks like a
            % measurement; NaN says 'not measurable' instead.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',256, 'SamplingFrequency',256);
            lam = j.Lambda;  om = j.Omega;
            C = exp(-(om - 1e4*sqrt(lam)).^2 / (2*(0.05*max(om))^2));
            d = tc.verifyWarningFree(@() dispersion(j, C));
            tc.verifyTrue(isnan(d.slope));
        end

    end
end

% Author: Diellor Basha, 2026
