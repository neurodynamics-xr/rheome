classdef tGfbTransform < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function seed(~), rng(0); end
    end

    methods (Test)

        function wtMatchesFrameAnalysis(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            b  = struct('Phi',fx.Phi, 'Lambda',lam, 'Mass',fx.Mass, 'nV',fx.nV);
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',7, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            F = randn(fx.nV, 5);
            tc.verifyEqual(wt(g, F), rheome.filters.frame_analysis(b, F, f), 'AbsTol', 1e-10);
        end

        function iwtMatchesFrameSynthesis(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            b  = struct('Phi',fx.Phi, 'Lambda',lam, 'Mass',fx.Mass, 'nV',fx.nV);
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',7, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            F = randn(fx.nV, 3);
            W = wt(g, F);
            tc.verifyEqual(iwt(g, W, 'tight'), ...
                rheome.filters.frame_synthesis(b, W, f, 'tight'), 'AbsTol', 1e-10);
            tc.verifyEqual(iwt(g, W, 'dual'), ...
                rheome.filters.frame_synthesis(b, W, f, 'dual'),  'AbsTol', 1e-10);
        end

        function dualReconstructsExactlyForEveryFamily(tc)
            % THE headline property: 'dual' inverts ANY bank, not just a tight one.
            fx = gfbFixture();  lam = fx.Lambda;
            T  = rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam);
            F  = fx.Phi * randn(numel(lam), 3);        % in-basis, so nothing is lost
            for fam = {'mexhat','itersine','heat'}
                g = rheome.graphfilterbank(lam, 'Wavelet', fam{1}, 'Transform', T);
                tc.verifyEqual(iwt(g, wt(g, F), 'dual'), F, 'AbsTol', 1e-8, ...
                    sprintf('dual reconstruction, %s', fam{1}));
            end
        end

        function tightReconstructsUpToAForItersine(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Wavelet','itersine', 'NumFilters',10, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            F = fx.Phi * randn(numel(lam), 2);
            A = framebounds(g, lam).A;
            tc.verifyEqual(iwt(g, wt(g, F), 'tight'), A*F, 'RelTol', 1e-6);
        end

        function trailingDimensionsPassThrough(tc)
            % Time is the SIGNAL's axis, not the transform's.
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            tc.verifySize(wt(g, randn(fx.nV, 1)),  [fx.nV, 1,  g.NumMembers]);
            tc.verifySize(wt(g, randn(fx.nV, 37)), [fx.nV, 37, g.NumMembers]);
        end

        function wtWithoutTransformErrorsClearly(tc)
            g = rheome.graphfilterbank(100);
            tc.verifyError(@() wt(g, randn(10,2)), 'graphfilterbank:noTransform');
        end

        function iwtRejectsAnUnknownMode(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            W  = wt(g, randn(fx.nV, 2));
            tc.verifyError(@() iwt(g, W, 'nonsense'), 'graphfilterbank:mode');
        end

    end
end

% Author: Diellor Basha, 2026
