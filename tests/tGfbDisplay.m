classdef tGfbDisplay < matlab.unittest.TestCase

    methods (Test)

        function spectralResponseReturnsGainsAndWavenumbers(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            [H, k] = spectralResponse(g);
            tc.verifySize(H, [numel(k), g.NumMembers]);
            tc.verifyTrue(all(k >= 0));
            tc.verifyEqual(max(k), sqrt(100), 'RelTol', 1e-12);
        end

        function spectralResponsePlotsWhenNoOutput(tc)
            g = rheome.graphfilterbank(100);
            f = figure('Visible','off');
            c = onCleanup(@() close(f));
            spectralResponse(g);
            tc.verifyNotEmpty(get(gca, 'Children'));
        end

        function dispMentionsBoundsAndUnusableMembers(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            s = evalc('disp(g)');
            tc.verifySubstring(s, 'mexhat');
            tc.verifySubstring(s, 'A =');
            tc.verifySubstring(s, 'unusable');
        end

        function dispSaysWhenNoTransformAttached(tc)
            s = evalc('disp(rheome.graphfilterbank(100))');
            tc.verifySubstring(s, 'no Transform');
        end

        function dispIsQuietAboutTransformWhenAttached(tc)
            fx = gfbFixture();
            g  = rheome.graphfilterbank(fx.Lambda, ...
                     'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, fx.Lambda));
            s = evalc('disp(g)');
            tc.verifyEmpty(strfind(s, 'no Transform'));
        end

    end
end

% Author: Diellor Basha, 2026
