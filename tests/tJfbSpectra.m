classdef tJfbSpectra < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function seed(~), rng(0); end
    end

    methods (Test)

        function scalogramIsMemberByFrequency(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13]);
            C  = randn(numel(j.Lambda), j.NumOmega);
            tc.verifySize(scalogram(j, C), [j.NumMembers, j.NumOmega]);
        end

        function scalogramEqualsTheDirectPerMemberEnergy(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13]);
            C  = randn(numel(j.Lambda), j.NumOmega);
            E  = scalogram(j, C);
            for m = 1:j.NumMembers
                direct = sum(abs(C .* jointfilters(j,m)).^2, 1);
                tc.verifyEqual(E(m,:), direct, 'RelTol', 1e-10, sprintf('member %d', m));
            end
        end

        function scalogramTotalMatchesTheCoefficientTotal(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            C  = randn(numel(j.Lambda), j.NumOmega);
            Y  = wt(j, C);
            tc.verifyEqual(sum(scalogram(j,C), 'all'), sum(abs(Y).^2, 'all'), 'RelTol', 1e-10);
        end

        function scalogramIsNonNegativeAndReal(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            C  = randn(numel(j.Lambda), j.NumOmega) + 1i*randn(numel(j.Lambda), j.NumOmega);
            E  = scalogram(j, C);
            tc.verifyTrue(isreal(E));
            tc.verifyTrue(all(E(:) >= 0));
        end

        function scalogramWorksBelowTheGuard(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'MaxBytes', 1000);
            C  = randn(numel(j.Lambda), j.NumOmega);
            tc.verifySize(scalogram(j, C), [j.NumMembers, j.NumOmega]);
        end

    end
end

% Author: Diellor Basha, 2026
