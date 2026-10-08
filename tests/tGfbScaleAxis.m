classdef tGfbScaleAxis < matlab.unittest.TestCase

    methods (Test)

        function centersMatchFiltersFrame(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters', 7);
            tc.verifyEqual(centerWavenumbers(g), f.Centers, 'RelTol', 1e-10);
        end

        function scalesAndWidthsMatchFiltersFrame(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters', 7);
            tc.verifyEqual(scales(g), f.t,     'RelTol', 1e-12);
            tc.verifyEqual(widths(g), f.Sigma, 'RelTol', 1e-12);
        end

        function widthIsExactlySqrt2t(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            tc.verifyEqual(widths(g), sqrt(2*scales(g)), 'RelTol', 1e-14);
        end

        function wavelengthIsTwoPiOverWavenumber(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            tc.verifyEqual(wavelengths(g), 2*pi./centerWavenumbers(g), 'RelTol', 1e-12);
        end

        function massLostMatchesFiltersFrame(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters', 7);
            tc.verifyEqual(g.MassLost, f.MassLost, 'RelTol', 1e-8);
            tc.verifyEqual(g.Usable,   f.Usable);
        end

        function finestScaleIsTheUsableFloor(tc)
            g = rheome.graphfilterbank(400);
            tc.verifyEqual(g.FinestScale, sqrt(2*4.744/400), 'RelTol', 1e-12);
        end

        function powerbwReturnsWavenumberBands(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            bw = powerbw(g);
            tc.verifySize(bw, [g.NumMembers, 2]);
            tc.verifyTrue(all(bw(:,2) >= bw(:,1)));
        end

        function qfactorIsDimensionless(tc)
            % Q is a ratio, so scaling the whole lambda axis must leave it unchanged.
            gA = rheome.graphfilterbank(100,   'Wavelet','mexhat', 'VoicesPerOctave',3);
            gB = rheome.graphfilterbank(10000, 'Wavelet','mexhat', 'VoicesPerOctave',3);
            tc.verifyEqual(qfactor(gB), qfactor(gA), 'RelTol', 1e-6);
        end

    end
end

% Author: Diellor Basha, 2026
