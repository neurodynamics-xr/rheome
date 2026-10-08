classdef tGfbBounds < matlab.unittest.TestCase

    methods (Test)

        function matchesFiltersFrameBounds(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet', 'mexhat', 'NumFilters', 7);
            ref = rheome.filters.frame_bounds(f, lam);
            b   = framebounds(g, lam);
            tc.verifyEqual(b.A, ref.A, 'AbsTol', 0);
            tc.verifyEqual(b.B, ref.B, 'AbsTol', 0);
            tc.verifyEqual(b.Tightness, ref.Tightness, 'AbsTol', 0);
        end

        function itersineIsTight(tc)
            g = rheome.graphfilterbank(100, 'Wavelet', 'itersine', 'NumFilters', 10);
            b = framebounds(g);
            tc.verifyEqual(b.Tightness, 1, 'RelTol', 1e-6);
            tc.verifyTrue(isframetight(g));
        end

        function mexhatIsNotTight(tc)
            g = rheome.graphfilterbank(100, 'Wavelet', 'mexhat', 'VoicesPerOctave', 3);
            tc.verifyFalse(isframetight(g));
        end

        function boundsFollowTheSuppliedSpectrum(tc)
            % Bounds are a property of the GRID the transform will see, not of the
            % continuum -- so two spectra on one bank must give different bounds.
            g  = rheome.graphfilterbank(100, 'Wavelet','mexhat', 'VoicesPerOctave',1);
            Ai = framebounds(g).A;                              % interval
            A1 = framebounds(g, linspace(1, 100, 50)').A;
            A2 = framebounds(g, linspace(40, 60, 50)').A;       % inside one member's band
            tc.verifyNotEqual(A1, A2);
            % A spectrum is a SUBSET of the interval, so it can never see a worse minimum.
            tc.verifyGreaterThanOrEqual(A1, Ai - 1e-12);
            tc.verifyGreaterThanOrEqual(A2, Ai - 1e-12);
        end

        function spectrumBackedBankIsAlreadyExact(tc)
            % Retaining Lambda is the whole reason a vector may be passed: the no-arg
            % call then needs no second argument to be exact.
            fx = gfbFixture();
            gs = rheome.graphfilterbank(fx.Lambda, 'Wavelet', 'mexhat');
            tc.verifyEqual(framebounds(gs).A, framebounds(gs, fx.Lambda).A, 'AbsTol', 0);
        end

        function degenerateFrameWarns(tc)
            % Family-aware, driven by computed bounds -- not a hardcoded 1.5/octave.
            % Coarse-only limits at a high lambda_max: every member peaks near
            % lambda ~ 10 while the axis runs to 1e4, so 95% of it is uncovered.
            tc.verifyWarning( ...
                @() rheome.graphfilterbank(1e4, 'Wavelet','mexhat', 'SizeLimits',[0.3 0.5]), ...
                'graphfilterbank:degenerateFrame');
        end

        function healthyFrameDoesNotWarn(tc)
            % Including heat, whose B/A ~ 56 would trip a naive tightness threshold.
            tc.verifyWarningFree(@() rheome.graphfilterbank(100, 'Wavelet','itersine', 'NumFilters',10));
            tc.verifyWarningFree(@() rheome.graphfilterbank(100, 'Wavelet','mexhat'));
            tc.verifyWarningFree(@() rheome.graphfilterbank(100, 'Wavelet','heat'));
        end

    end
end

% Author: Diellor Basha, 2026
