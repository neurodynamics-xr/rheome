classdef tJfbFrame < matlab.unittest.TestCase

    methods (Test)

        function frameOperatorIsTheSumOfSquaredMagnitudes(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13]);
            b  = framebounds(j);
            S  = zeros(numel(j.Lambda), j.NumOmega);
            for m = 1:j.NumMembers, S = S + abs(jointfilters(j,m)).^2; end
            tc.verifyEqual(b.S, S, 'AbsTol', 1e-12);
            tc.verifyEqual(b.A, min(S(:)), 'AbsTol', 0);
            tc.verifyEqual(b.B, max(S(:)), 'AbsTol', 0);
        end

        function boundsComeWithCoverageAndWorstPoint(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            b  = framebounds(j);
            tc.verifyTrue(isfield(b,'Uncovered') && isfield(b,'Worst'));
            tc.verifySize(b.Worst, [1 2]);
        end

        function anAllPassBankOverATightGraphBankIsTight(tc)
            % itersine tiles the lambda axis tightly and the default time member is an
            % all-pass, so the joint frame operator is flat.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            tc.verifyEqual(framebounds(j).Tightness, 1, 'RelTol', 1e-6);
            tc.verifyTrue(isframetight(j));
        end

        function aNarrowBandLeavesThePlaneUncovered(tc)
            % Outside the band every member is zero, so A collapses.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Bands', [8 13]);
            b = framebounds(j);
            tc.verifyLessThan(b.A, 1e-12);
            tc.verifyGreaterThan(b.Uncovered, 0);
            tc.verifyFalse(isframetight(j));
        end

        function boundsNeverMaterialiseTheBank(tc)
            % framebounds must work below the guard -- this is why it loops.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'MaxBytes', 1000);
            tc.verifyNotEmpty(framebounds(j));
        end

    end
end

% Author: Diellor Basha, 2026
