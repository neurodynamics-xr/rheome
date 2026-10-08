classdef tGfbFilters < matlab.unittest.TestCase

    methods (Test)

        function gainsOnSuppliedSpectrum(tc)
            g   = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            lam = linspace(0, 100, 257)';
            H   = graphfilters(g, 'Lambda', lam);
            tc.verifySize(H, [257, g.NumMembers]);
            tc.verifyTrue(all(isfinite(H(:))));
        end

        function gainsDefaultToRetainedSpectrum(tc)
            fx = gfbFixture();
            g  = rheome.graphfilterbank(fx.Lambda);
            tc.verifySize(graphfilters(g), [numel(fx.Lambda), g.NumMembers]);
        end

        function gainsFallBackToDenseGridWithoutSpectrum(tc)
            g = rheome.graphfilterbank(100);
            H = graphfilters(g);
            tc.verifyEqual(size(H,2), g.NumMembers);
            tc.verifyGreaterThanOrEqual(size(H,1), 256);
        end

        function dualReconstructsUnityWhereCovered(tc)
            % sum_m gd_m * g_m = 1 wherever S > 0.
            g   = rheome.graphfilterbank(100, 'Wavelet','mexhat', 'VoicesPerOctave', 4);
            lam = linspace(0.1, 100, 512)';
            H   = graphfilters(g, 'Lambda', lam);
            Hd  = graphfilters(g, 'Lambda', lam, 'Dual', true);
            S   = sum(H.^2, 2);
            ok  = S > 1e-12*max(S);
            tc.verifyEqual(sum(Hd(ok,:).*H(ok,:), 2), ones(nnz(ok),1), 'AbsTol', 1e-12);
        end

        function dualIsZeroWhereUncovered(tc)
            g   = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            lam = [0; 1e6];                       % 1e6 is far outside the bank
            Hd  = graphfilters(g, 'Lambda', lam, 'Dual', true);
            tc.verifyEqual(Hd(2,:), zeros(1, g.NumMembers), 'AbsTol', 0);
        end

        function gainReturnsOneMemberHandle(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            h = gain(g, 2);
            tc.verifyClass(h, 'function_handle');
            H = graphfilters(g, 'Lambda', [1;2;3]);
            tc.verifyEqual(h([1;2;3]), H(:,2), 'AbsTol', 0);
            tc.verifyError(@() gain(g, 0), 'graphfilterbank:member');
        end

        function unknownTypeIsRefused(tc)
            g = rheome.graphfilterbank(100);
            tc.verifyError(@() graphfilters(g, 'Type', 'nonsense'), 'graphfilterbank:type');
        end

    end
end

% Author: Diellor Basha, 2026
