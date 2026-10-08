classdef tJfbTransform < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function seed(~), rng(0); end
    end

    methods (Test)

        function wtIsAnElementwiseMultiplyPerMember(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13]);
            C  = randn(numel(j.Lambda), j.NumOmega) + 1i*randn(numel(j.Lambda), j.NumOmega);
            Y  = wt(j, C);
            tc.verifySize(Y, [numel(j.Lambda), j.NumOmega, j.NumMembers]);
            for m = 1:j.NumMembers
                tc.verifyEqual(Y(:,:,m), C .* jointfilters(j,m), 'AbsTol', 1e-12);
            end
        end

        function dualReconstructsExactlyWhereCovered(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            C  = randn(numel(j.Lambda), j.NumOmega) + 1i*randn(numel(j.Lambda), j.NumOmega);
            tc.verifyEqual(iwt(j, wt(j, C), 'dual'), C, 'AbsTol', 1e-8);
        end

        function dualReconstructsANonSeparableBankToo(tc)
            % Speed-selective members are not rank-1; the dual must still invert.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                     'JointKernels', {@(l,w) 0.5 + 0.5*exp(-(w - 5*sqrt(l)).^2/(2*20^2))});
            C = randn(numel(j.Lambda), j.NumOmega);
            tc.verifyEqual(iwt(j, wt(j, C), 'dual'), C, 'AbsTol', 1e-8);
        end

        function tightGivesSTimesC(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            C  = randn(numel(j.Lambda), j.NumOmega);
            S  = framebounds(j).S;
            tc.verifyEqual(iwt(j, wt(j, C), 'tight'), S .* C, 'AbsTol', 1e-10);
        end

        function uncoveredContentIsLostNotAmplified(tc)
            % Outside a band the dual is zero, so iwt returns zero there -- it must not
            % divide by a vanishing S and blow up.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Bands',[8 13]);
            C  = ones(numel(j.Lambda), j.NumOmega);
            R  = iwt(j, wt(j, C), 'dual');
            tc.verifyTrue(all(isfinite(R(:))));
            out = j.Frequencies < 8 | j.Frequencies > 13;
            tc.verifyEqual(R(:, out), zeros(numel(j.Lambda), nnz(out)), 'AbsTol', 1e-12);
        end

        function iwtRejectsAnUnknownMode(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            Y  = wt(j, randn(numel(j.Lambda), j.NumOmega));
            tc.verifyError(@() iwt(j, Y, 'nonsense'), 'jointfilterbank:mode');
        end

        function wtChecksTheSpectrumSize(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            tc.verifyError(@() wt(j, randn(3,4)), 'jointfilterbank:size');
        end

    end
end

% Author: Diellor Basha, 2026
