classdef tGfbChebyshev < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function seed(~), rng(0); end
    end

    methods (Test)

        function chebyshevApproximatesTheEigenTransform(tc)
            fx = gfbFixture();  lam = fx.Lambda;  lmaxTrunc = max(lam);
            Te = rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam);
            Tc = rheome.graphtransform.chebyshev(fx.L, [], 'Order', 160, 'Mass', fx.M);
            % Same DESIGN on both, so only the application route differs.
            ge = rheome.graphfilterbank(lmaxTrunc, 'Wavelet','mexhat','NumFilters',5,'Transform',Te);
            gc = rheome.graphfilterbank(lmaxTrunc, 'Wavelet','mexhat','NumFilters',5,'Transform',Tc);
            % Band-limited field: it lies wholly in the retained span, so the truncated
            % eigenbasis is exact there and any difference is Chebyshev error alone.
            F  = fx.Phi * randn(numel(lam), 1);
            We = wt(ge, F);  Wc = wt(gc, F);
            tc.verifyLessThan(norm(Wc(:)-We(:))/norm(We(:)), 5e-2);
        end

        function higherOrderIsMoreAccurate(tc)
            fx = gfbFixture();  lam = fx.Lambda;  lmaxTrunc = max(lam);
            Te = rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam);
            ge = rheome.graphfilterbank(lmaxTrunc,'Wavelet','mexhat','NumFilters',5,'Transform',Te);
            F  = fx.Phi * randn(numel(lam), 1);
            ref = wt(ge, F);
            ords = [40 160];  err = zeros(1,2);
            for i = 1:2
                Tc = rheome.graphtransform.chebyshev(fx.L, [], 'Order', ords(i), 'Mass', fx.M);
                gc = rheome.graphfilterbank(lmaxTrunc,'Wavelet','mexhat','NumFilters',5,'Transform',Tc);
                Wc = wt(gc, F);
                err(i) = norm(Wc(:) - ref(:)) / norm(ref(:));
            end
            tc.verifyLessThan(err(2), err(1));
        end

        function fromOperatorNeedsNoEigenvectors(tc)
            % The whole point: lmax alone, obtained by Lanczos.
            fx = gfbFixture();
            g  = rheome.graphfilterbank.fromOperator(fx.L, 'Mass', fx.M, 'Wavelet','mexhat');
            tc.verifyGreaterThan(g.SpectralRange, 0);
            tc.verifyNotEmpty(g.Transform);
            tc.verifySize(wt(g, randn(fx.nV, 2)), [fx.nV, 2, g.NumMembers]);
        end

        function fromOperatorForwardsConstructorOptions(tc)
            fx = gfbFixture();
            g  = rheome.graphfilterbank.fromOperator(fx.L, 'Mass', fx.M, ...
                     'Wavelet','itersine', 'NumFilters', 6);
            tc.verifyEqual(g.Wavelet, 'itersine');
            tc.verifyEqual(g.NumMembers, 6);
        end

        function truncatedLmaxIsCaughtNotSilentlyDivergent(tc)
            % Reusing a truncated basis's lambda_max makes the recursion diverge to
            % 1e31 and then NaN. It must be refused loudly, not returned.
            fx = gfbFixture();
            tc.verifyWarning( ...
                @() rheome.graphtransform.chebyshev(fx.L, max(fx.Lambda), 'Mass', fx.M), ...
                'graphtransform:chebyshev:lmax');
        end

        function estimatedLmaxDoesNotWarn(tc)
            fx = gfbFixture();
            tc.verifyWarningFree(@() rheome.graphtransform.chebyshev(fx.L, [], 'Mass', fx.M));
        end

        function chebyshevWaveletSpectraMatchEigenSpectra(tc)
            % vertexSpectrum must take the .filter route and agree on the WAVELET
            % members. Member 1 is excluded deliberately -- see the next test.
            [Ee, Ec] = i_spectra(160);
            tc.verifyEqual(Ec(2:end,:), Ee(2:end,:), 'RelTol', 0.05);
        end

        function lowPassIsTheChebyshevWeakSpotAndConverges(tc)
            % The scaling function is a QUARTIC roll-off at lambda ~ 186 while the
            % polynomial is fitted over [0, 4.7e5] -- a feature spanning 0.04% of the
            % interval, which a polynomial resolves only at very high order. Measured:
            % 35% error at order 160, 5.6% at order 320. That is convergence, not a
            % defect, but it means the eigen transform is the right route when the
            % low-pass band matters.
            [Ee1, Ec1] = i_spectra(160);
            [Ee2, Ec2] = i_spectra(320);
            e160 = abs(Ec1(1,1) - Ee1(1,1)) / Ee1(1,1);
            e320 = abs(Ec2(1,1) - Ee2(1,1)) / Ee2(1,1);
            tc.verifyGreaterThan(e160, 0.1);      % genuinely inaccurate at 160
            tc.verifyLessThan(e320, e160 / 2);    % and genuinely converging
        end

    end
end

function [Ee, Ec] = i_spectra(ord)
    rng(0);
    fx = gfbFixture();  lam = fx.Lambda;  lt = max(lam);
    Te = rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam);
    Tc = rheome.graphtransform.chebyshev(fx.L, [], 'Order', ord, 'Mass', fx.M);
    ge = rheome.graphfilterbank(lt, 'Wavelet','mexhat', 'NumFilters',5, 'Transform',Te);
    gc = rheome.graphfilterbank(lt, 'Wavelet','mexhat', 'NumFilters',5, 'Transform',Tc);
    F  = fx.Phi * randn(numel(lam), 2);
    Ee = vertexSpectrum(ge, F);
    Ec = vertexSpectrum(gc, F);
end

% Author: Diellor Basha, 2026
