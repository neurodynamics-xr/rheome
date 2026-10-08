classdef tJfbAxes < matlab.unittest.TestCase
% The temporal grid, the half/boundary conventions, and compatibility checking.

    methods (Test)

        function signalLengthDerivesThePositiveHalf(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength', 64, 'SamplingFrequency', 64);
            tc.verifyEqual(j.NumOmega, 33);              % bins 0..32 inclusive
            tc.verifyEqual(j.Frequencies(1),   0,  'AbsTol', 0);
            tc.verifyEqual(j.Frequencies(end), 32, 'RelTol', 1e-12);
            tc.verifyEqual(j.Df, 1, 'RelTol', 1e-12);
        end

        function fullHalfKeepsEveryBin(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Half','full');
            tc.verifyEqual(j.NumOmega, 64);
        end

        function samplingFrequencyIsPureLabelling(tc)
            fx = jfbFixture();
            j1 = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',1);
            j2 = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64);
            tc.verifyEqual(j2.Frequencies, 64*j1.Frequencies, 'RelTol', 1e-12);
            tc.verifyEqual(j2.NumOmega, j1.NumOmega);
        end

        function explicitFrequenciesAreHonoured(tc)
            fx = jfbFixture();
            f  = linspace(8, 13, 17);
            j  = rheome.jointfilterbank(fx.gfb, 'Frequencies', f);
            tc.verifyEqual(j.Frequencies, f, 'RelTol', 1e-12);
            tc.verifyEqual(j.NumOmega, 17);
        end

        function omegaIsTwoPiF(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64);
            tc.verifyEqual(j.Omega, 2*pi*j.Frequencies, 'RelTol', 1e-12);
        end

        function temporalGridIsRequired(tc)
            fx = jfbFixture();
            tc.verifyError(@() rheome.jointfilterbank(fx.gfb), 'jointfilterbank:temporalGrid');
        end

        function bothGridFormsAtOnceIsRefused(tc)
            fx = jfbFixture();
            tc.verifyError(@() rheome.jointfilterbank(fx.gfb, 'SignalLength',64, ...
                'Frequencies', 1:10), 'jointfilterbank:temporalGrid');
        end

        function graphSideMayBeALambdaMax(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(max(fx.Lambda), 'SignalLength', 64);
            tc.verifyClass(j.GraphBank, 'rheome.graphfilterbank');
        end

        function axesCarriesBothAxesAndTheConventions(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64);
            ax = axes(j);
            tc.verifyEqual(ax.nOmega, 33);
            tc.verifyEqual(ax.half, 'positive');
            tc.verifyEqual(ax.boundary, 'ring');
            tc.verifyEqual(ax.k, sqrt(ax.lambda), 'RelTol', 1e-12);
        end

        function iscompatibleCatchesWhatASizeCheckCannot(tc)
            % Two bands can retain the SAME bin count and still not belong together.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'Frequencies', linspace(8,13,17));
            other = rheome.jointfilterbank(fx.gfb, 'Frequencies', linspace(20,30,17));
            [ok, why] = iscompatible(j, axes(other));
            tc.verifyFalse(ok);
            tc.verifyNotEmpty(why);
        end

        function iscompatibleCatchesHalfAndBoundary(tc)
            fx = jfbFixture();
            a = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Half','positive');
            b = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Half','full');
            c = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Boundary','path');
            tc.verifyFalse(iscompatible(a, axes(b)));
            tc.verifyFalse(iscompatible(a, axes(c)));
            tc.verifyTrue(iscompatible(a, axes(a)));
        end

    end
end

% Author: Diellor Basha, 2026
