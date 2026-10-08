classdef tGfbSpectra < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function seed(~), rng(0); end
    end

    methods (Test)

        function vertexSpectrumMatchesFrameScalogram(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            b  = struct('Phi',fx.Phi, 'Lambda',lam, 'Mass',fx.Mass, 'nV',fx.nV);
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',7, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            F   = randn(fx.nV, 6);
            ref = rheome.filters.frame_scalogram(b, F, f);
            tc.verifyEqual(vertexSpectrum(g, F), ref.energy, 'RelTol', 1e-8);
        end

        function vertexSpectrumShapeFollowsTheSignal(tc)
            % The [M x nT] scalogram falls out free, with no time axis in the class.
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            tc.verifySize(vertexSpectrum(g, randn(fx.nV, 1)),  [g.NumMembers, 1]);
            tc.verifySize(vertexSpectrum(g, randn(fx.nV, 40)), [g.NumMembers, 40]);
        end

        function vertexSpectrumEqualsDirectBandEnergy(tc)
            % Parseval: coefficient-space energy IS the field energy.
            fx = gfbFixture();  lam = fx.Lambda;
            T  = rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',5, 'Transform',T);
            F  = fx.Phi * randn(numel(lam), 2);
            W  = wt(g, F);
            E  = vertexSpectrum(g, F);
            for m = 1:g.NumMembers
                direct = sum(abs(W(:,:,m)).^2, 1);
                tc.verifyEqual(E(m,:), direct, 'RelTol', 1e-6, sprintf('member %d', m));
            end
        end

        function scaleSpectrumIsAPerVertexMap(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            F  = randn(fx.nV, 4);
            S  = scaleSpectrum(g, F);
            tc.verifySize(S, [fx.nV, 4]);
            tc.verifyTrue(all(S(:) >= 0));
        end

        function bothMarginalsShareOneTotal(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',5, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            F = fx.Phi * randn(numel(lam), 3);
            tc.verifyEqual(sum(vertexSpectrum(g,F), 1), sum(scaleSpectrum(g,F), 1), ...
                'RelTol', 1e-6);
        end

    end
end

% Author: Diellor Basha, 2026
