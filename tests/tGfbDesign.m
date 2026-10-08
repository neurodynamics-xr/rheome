classdef tGfbDesign < matlab.unittest.TestCase
% Bank design, and bit-for-bit regression against the validated rheome.filters.frame.

    methods (Test)

        function mexhatMatchesFiltersFrame(tc)
            fx = gfbFixture();
            lam = fx.Lambda;
            for Nf = [4 7 12]
                f = rheome.filters.frame('mexhat', Nf, lam, 'warn', false);
                g = rheome.graphfilterbank(lam, 'Wavelet', 'mexhat', 'NumFilters', Nf);
                tc.verifyEqual(g.NumMembers, numel(f.g), ...
                    sprintf('member count, Nf = %d', Nf));
                Href = rheome.filters.frame_gains(f, lam);
                Hnew = graphfilters(g);
                tc.verifyEqual(Hnew, Href, 'AbsTol', 0, ...
                    sprintf('gains must match EXACTLY, Nf = %d', Nf));
            end
        end

        function itersineMatchesFiltersFrame(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            f = rheome.filters.frame('itersine', 8, lam, 'warn', false);
            g = rheome.graphfilterbank(lam, 'Wavelet', 'itersine', 'NumFilters', 8);
            tc.verifyEqual(graphfilters(g), rheome.filters.frame_gains(f, lam), 'AbsTol', 0);
        end

        function heatMatchesFiltersFrame(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            f = rheome.filters.frame('heat', 6, lam, 'warn', false);
            g = rheome.graphfilterbank(lam, 'Wavelet', 'heat', 'NumFilters', 6);
            tc.verifyEqual(graphfilters(g), rheome.filters.frame_gains(f, lam), 'AbsTol', 0);
        end

        function derivedCountMatchesFiltersFrame(tc)
            % Nf = [] -> derived at VoicesPerOctave members/octave, as frame.m does.
            fx = gfbFixture();  lam = fx.Lambda;
            f = rheome.filters.frame('mexhat', [], lam, 'perOctave', 3, 'warn', false);
            g = rheome.graphfilterbank(lam, 'Wavelet', 'mexhat', 'VoicesPerOctave', 3);
            tc.verifyEqual(g.NumMembers, numel(f.g));
            tc.verifyEqual(graphfilters(g), rheome.filters.frame_gains(f, lam), 'AbsTol', 0);
        end

        function voicesPerOctaveChangesDensityNotReach(tc)
            % Voices controls DENSITY only; the scale range must not move.
            g3 = rheome.graphfilterbank(100, 'VoicesPerOctave', 3);
            g8 = rheome.graphfilterbank(100, 'VoicesPerOctave', 8);
            tc.verifyEqual(g8.ScaleLimits, g3.ScaleLimits, 'RelTol', 1e-12);
            tc.verifyGreaterThan(g8.NumMembers, g3.NumMembers);
        end

        function voicesPerOctaveIsValidated(tc)
            tc.verifyError(@() rheome.graphfilterbank(100, 'VoicesPerOctave', 0),   'MATLAB:validators:mustBePositive');
            tc.verifyError(@() rheome.graphfilterbank(100, 'VoicesPerOctave', 49),  'graphfilterbank:voices');
            tc.verifyError(@() rheome.graphfilterbank(100, 'VoicesPerOctave', 2.5), 'MATLAB:validators:mustBeInteger');
        end

        function explicitNumFiltersReportsAchievedVoices(tc)
            g = rheome.graphfilterbank(100, 'Wavelet', 'mexhat', 'NumFilters', 5);
            tc.verifyEqual(g.NumMembers, 6);            % 5 wavelets + 1 scaling function
            tc.verifyTrue(isfinite(g.AchievedVoicesPerOctave));
        end

        function sizeLimitsMakeBankInvariantToLambdaMax(tc)
            % The K-invariance property, asserted explicitly.
            s = [6e-3 60e-3];
            gA = rheome.graphfilterbank(400,  'SizeLimits', s, 'VoicesPerOctave', 3);
            gB = rheome.graphfilterbank(1600, 'SizeLimits', s, 'VoicesPerOctave', 3);
            tc.verifyEqual(gB.ScaleLimits, gA.ScaleLimits, 'RelTol', 1e-12);
            tc.verifyEqual(scales(gB), scales(gA), 'RelTol', 1e-12);
        end

        function lpFactorDefaultMovesCoarseEnd(tc)
            % The documented DEFECT of the default. Asserted so it cannot regress silently.
            gA = rheome.graphfilterbank(400);
            gB = rheome.graphfilterbank(1600);
            tc.verifyLessThan(gB.SizeLimits(2), gA.SizeLimits(2));
        end

    end
end

% Author: Diellor Basha, 2026
