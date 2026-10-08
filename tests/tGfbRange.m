classdef tGfbRange < matlab.unittest.TestCase
% Spectral range resolution: the three limit views, the LPFactor default, signed refusal.

    methods (Test)

        function scalarRangeUsesLPFactor(tc)
            g = rheome.graphfilterbank(100);
            tc.verifyEqual(g.SpectralRange, 100);
            % lminEff = lmax/LPFactor = 5;  tmax = 2/lminEff;  tmin = 1/lmax ('basis')
            tc.verifyEqual(g.ScaleLimits, [1/100, 2/5], 'RelTol', 1e-12);
        end

        function explicitRangeIsHonoured(tc)
            g = rheome.graphfilterbank([2 100]);
            tc.verifyEqual(g.ScaleLimits, [1/100, 2/2], 'RelTol', 1e-12);
        end

        function spectrumVectorRetainsExactLambda(tc)
            fx = gfbFixture();
            g  = rheome.graphfilterbank(fx.Lambda);
            tc.verifyEqual(g.SpectralRange, max(fx.Lambda), 'RelTol', 1e-12);
            tc.verifyTrue(g.HasSpectrum);
        end

        function threeLimitViewsAgree(tc)
            % Setting any one view must produce an identical bank. sigma = sqrt(2t).
            t  = [2e-5 8e-4];
            gT = rheome.graphfilterbank(100, 'ScaleLimits', t);
            gS = rheome.graphfilterbank(100, 'SizeLimits',  sqrt(2*t));
            gL = rheome.graphfilterbank(100, 'SpectralLimits', [2/t(2), 1/t(1)]);
            tc.verifyEqual(gS.ScaleLimits, gT.ScaleLimits, 'RelTol', 1e-12);
            tc.verifyEqual(gL.ScaleLimits, gT.ScaleLimits, 'RelTol', 1e-12);
        end

        function sizeLimitsRoundTrip(tc)
            s = [6e-3 60e-3];
            g = rheome.graphfilterbank(100, 'SizeLimits', s);
            tc.verifyEqual(g.SizeLimits, s, 'RelTol', 1e-12);
        end

        function fineLimitUsableMovesTminOnly(tc)
            gb = rheome.graphfilterbank(100, 'FineLimit', 'basis');
            gu = rheome.graphfilterbank(100, 'FineLimit', 'usable');
            tc.verifyEqual(gu.ScaleLimits(1), 4.744/100, 'RelTol', 1e-12);
            tc.verifyEqual(gu.ScaleLimits(2), gb.ScaleLimits(2), 'RelTol', 1e-12);
        end

        function zeroLambdaMinFallsBackToLPFactor(tc)
            % lambda_min = 0 is the DC mode and carries no scale information; using it
            % directly would give t_max = Inf.
            g0 = rheome.graphfilterbank([0 100]);
            gs = rheome.graphfilterbank(100);
            tc.verifyEqual(g0.ScaleLimits, gs.ScaleLimits, 'RelTol', 1e-12);
            tc.verifyTrue(all(isfinite(g0.ScaleLimits)));
        end

        function signedSpectrumIsRefused(tc)
            % An unsquared Dirac has signed lambda. mexhat/heat diverge and
            % sqrt(lambda) is imaginary. Refuse loudly rather than return NaN scales.
            tc.verifyError(@() rheome.graphfilterbank([-3 100]), ...
                'graphfilterbank:signedSpectrum');
        end

        function degenerateRangeIsRefused(tc)
            tc.verifyError(@() rheome.graphfilterbank([100 2]), 'graphfilterbank:range');
            tc.verifyError(@() rheome.graphfilterbank(0),       'graphfilterbank:range');
        end

    end
end

% Author: Diellor Basha, 2026
