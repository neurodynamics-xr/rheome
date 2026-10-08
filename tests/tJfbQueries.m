classdef tJfbQueries < matlab.unittest.TestCase

    methods (Test)

        function centresAreOnePerMember(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Bands',[4 8; 8 13]);
            tc.verifySize(centerWavenumbers(j), [1, j.NumMembers]);
            tc.verifySize(centerFrequencies(j), [1, j.NumMembers]);
        end

        function centreFrequencyLandsInsideItsBand(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',256, 'SamplingFrequency',256, ...
                                 'Bands',[8 13]);
            cf = centerFrequencies(j);
            tc.verifyTrue(all(cf >= 8 & cf <= 13));
        end

        function centreFrequencyIsInHzAndScalesWithFs(tc)
            % Hz is the CORRECT word on the time axis -- unlike lambda, omega really is
            % a frequency. Doubling fs doubles the centre.
            fx = jfbFixture();
            j1 = rheome.jointfilterbank(fx.gfb, 'SignalLength',256, 'SamplingFrequency',128);
            j2 = rheome.jointfilterbank(fx.gfb, 'SignalLength',256, 'SamplingFrequency',256);
            tc.verifyEqual(centerFrequencies(j2), 2*centerFrequencies(j1), 'RelTol', 1e-10);
        end

        function centreWavenumberTracksTheGraphBank(tc)
            % Members sharing a graph factor share a centre wavenumber.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[4 8; 8 13]);
            cw = centerWavenumbers(j);
            for ig = 1:j.NumGraph
                sel = j.Index(:,1) == ig;
                tc.verifyEqual(max(cw(sel)) - min(cw(sel)), 0, 'AbsTol', 1e-9);
            end
        end

        function spectralResponseReturnsTheMemberOnTheTwoAxes(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64);
            [W, k, f] = spectralResponse(j, 1);
            tc.verifySize(W, [numel(k), numel(f)]);
            tc.verifyEqual(k, sqrt(j.Lambda), 'RelTol', 1e-12);
            tc.verifyEqual(f, j.Frequencies,  'RelTol', 1e-12);
        end

        function spectralResponsePlotsWhenNoOutput(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            fig = figure('Visible','off');
            c = onCleanup(@() close(fig));
            spectralResponse(j, 1);
            tc.verifyNotEmpty(get(gca, 'Children'));
        end

        function dispReportsShapeBoundsAndLaziness(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13]);
            s = evalc('disp(j)');
            tc.verifySubstring(s, 'jointfilterbank');
            tc.verifySubstring(s, 'members');
            tc.verifySubstring(s, 'A =');
            tc.verifySubstring(s, 'materialised');
        end

        function dispSaysWhenNoTransformAttached(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength', 64);
            s = evalc('disp(j)');
            tc.verifySubstring(s, 'no Transform');
        end

    end
end

% Author: Diellor Basha, 2026
