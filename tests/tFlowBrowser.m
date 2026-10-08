classdef tFlowBrowser < matlab.unittest.TestCase
% @flowbrowser -- construction and the programmatic API.
%
% The point of these is not to test graphics. It is that the callbacks ARE the app's
% programmatic API, so driving band/page/frame/mode from a script must work without
% synthesising UI events -- which is also how the app gets scripted and animated.
%
% Runs against a small synthetic anatomy (an icosphere) so it costs milliseconds rather
% than the ~40 s a real rheome.flowpage.prepare takes.
%
% Author: Diellor Basha, 2026

    properties
        B; App
        C = 8; Ks = 12; nT = 600; fs = 400;
    end

    methods (TestClassSetup)
        function makeBundle(tc)
            if ~usejava('desktop') && isempty(ver('matlab')), tc.assumeFail('no graphics'); end
            d = tempname; mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            f = fullfile(d, 'recording.mat');

            rng(5);
            F = randn(tc.C, tc.nT);                                   %#ok<NASGU>
            Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;                 %#ok<NASGU>
            ChannelFlag = ones(tc.C,1);                               %#ok<NASGU>
            ChannelName = compose('C%d', 1:tc.C);                     %#ok<NASGU>
            ChannelType = repmat({'MEG'}, 1, tc.C);                   %#ok<NASGU>
            nCh = tc.C; nT = tc.nT; Comment = 's'; source = '';       %#ok<NASGU>
            save(f, 'F','Time','sfreq','ChannelFlag','ChannelName','ChannelType', ...
                 'nCh','nT','Comment','source', '-v7.3', '-nocompression');

            [Vx, Fc] = rheome.geom.icosphere(1);                             % 42 vertices
            nV = size(Vx,1);
            [Phi, ~] = qr(randn(nV, tc.Ks), 0);
            Lambda = linspace(1, 400, tc.Ks)';

            tc.B = struct('pager', rheome.pagedrecording(f,'PageLength',200,'Overlap',60,'Precision','single'), ...
                'Kc', randn(tc.Ks, tc.C), 'Phi', Phi, 'Lambda', Lambda, 'fs', tc.fs, ...
                'gfb', rheome.graphfilterbank(max(Lambda),'NumFilters',3), 'wVert', ones(nV,1), ...
                'surface', struct('Vertices',Vx,'Faces',Fc,'nV',nV), 'name','synthetic', ...
                'freqLimits',[10 120], 'voicesOct', 6, 'chSel', (1:tc.C).');
        end
    end

    methods (TestMethodTeardown)
        function shut(tc)
            if ~isempty(tc.App) && isvalid(tc.App), delete(tc.App); end
            tc.App = [];
        end
    end

    methods (Test)

        function opensWithTheExpectedTabs(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            tc.verifyTrue(isvalid(tc.App.Fig));
            tc.verifyEqual({tc.App.Tabs.Children.Title}, {'Cortex','Spectrum','Timeline'});
        end

        function landsOnACoreFrameNotInTheMargin(tc)
            % The margin belongs to the neighbouring pages and carries the wavelet cone;
            % opening on it would show contaminated frames as if they were the page's.
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 2);
            tc.verifyTrue(tc.App.Page.Core(tc.App.Frame));
        end

        function scrubbingMovesTheFrameWithoutRecomputingThePage(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            before = tc.App.Page.Coefficients;
            tc.App.onFrame(tc.App.Frame + 20);
            tc.verifyEqual(tc.App.Page.Coefficients, before);   % same page object contents
        end

        function frameIsClampedRatherThanErroring(tc)
            % A slider can overshoot; the app must clamp, not throw in a callback.
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            tc.App.onFrame(-5);
            tc.verifyEqual(tc.App.Frame, 1);
            tc.App.onFrame(1e6);
            tc.verifyEqual(tc.App.Frame, tc.App.Page.NumFrames);
        end

        function switchingBandRebuildsThePageAtItsOwnRate(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            rA = tc.App.Page.Rate;  nA = tc.App.Page.NumFilters;
            tc.App.onBand('gamma');
            tc.verifyEqual(tc.App.BandName, 'gamma');
            tc.verifyNotEqual(tc.App.Page.NumFilters, nA);
            tc.verifyGreaterThanOrEqual(tc.App.Page.Rate, rA);   % higher band, >= rate
        end

        function switchingModeDoesNotDisturbTheFrame(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            tc.App.onFrame(tc.App.Frame + 11);
            i = tc.App.Frame;
            tc.App.onMode('magnitude');
            tc.verifyEqual(tc.App.MapMode, 'magnitude');
            tc.verifyEqual(tc.App.Frame, i);
        end

        function stepButtonsAndKeyboardMoveTheFrameTogether(tc)
            % The step buttons, the slider, the keyboard and the play loop all go through
            % setFrame, so they cannot drift out of sync -- which is the failure that makes
            % a browser show one frame and label another.
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            i0 = tc.App.Frame;
            tc.App.stepFrame(1);
            tc.verifyEqual(tc.App.Frame, i0 + 1);
            tc.App.onKey(struct('Key','rightarrow','Modifier',{{}}));
            tc.verifyEqual(tc.App.Frame, i0 + 2);
            tc.App.onKey(struct('Key','rightarrow','Modifier',{{'shift'}}));
            tc.verifyEqual(tc.App.Frame, i0 + 12);          % shift = ten
            tc.App.onKey(struct('Key','leftarrow','Modifier',{{}}));
            tc.verifyEqual(tc.App.Frame, i0 + 11);
        end

        function theFrameLabelReportsTimeAndFlagsTheMargin(tc)
            % The Cortex tab has no time axis to shade, so without this the surfaces give no
            % indication that a frame is cone-contaminated.
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 2);
            tc.verifySubstring(tc.App.FrameLbl.Text, 's');
            tc.verifyFalse(contains(tc.App.FrameLbl.Text, 'margin'));   % opened on core
            tc.App.setFrame(1);                                          % into the margin
            tc.verifyFalse(tc.App.Page.Core(1));
            tc.verifySubstring(tc.App.FrameLbl.Text, 'margin');
        end

        function keyboardTogglesTheMapMode(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            tc.App.onKey(struct('Key','m','Modifier',{{}}));
            tc.verifyEqual(tc.App.MapMode, 'magnitude');
            tc.App.onKey(struct('Key','m','Modifier',{{}}));
            tc.verifyEqual(tc.App.MapMode, 'signed');
        end

        function theSensorStripTracksTheActiveFrame(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            tc.verifyTrue(isvalid(tc.App.SensorAx));
            tc.verifyEqual(numel(findobj(tc.App.SensorAx,'Type','line')), tc.App.NumSensorTraces);
            tc.App.setFrame(tc.App.Frame + 30);
            tc.verifyEqual(tc.App.SensorCursor.Value, tc.App.Page.Time(tc.App.Frame), 'AbsTol', 1e-12);
        end

        function theDivergingColormapHasBalancedArms(tc)
            % A naive blue-white-red gives the two arms unequal luminance, so equal-magnitude
            % CCW and CW rotation do not look equal and the eye reads a handedness bias that
            % is not in the data. Check the two ends are close in lightness.
            % ⚠ LUMINANCE MUST BE COMPUTED ON LINEARISED sRGB. Applying the Rec.709 weights
            % to gamma-ENCODED values is not luminance and reports these endpoints as
            % mismatched by 0.145 when in linear light they agree to 0.6%. Getting this
            % wrong would condemn a correctly balanced map.
            c = rheome.flowbrowser.coolwarm(256);
            tc.verifySize(c, [256 3]);
            tc.verifyTrue(all(c(:) >= 0 & c(:) <= 1));
            lin = c;
            hi  = c > 0.04045;
            lin(hi)  = ((c(hi) + 0.055) / 1.055) .^ 2.4;
            lin(~hi) = c(~hi) / 12.92;
            lum = lin * [0.2126; 0.7152; 0.0722];
            tc.verifyEqual(lum(1), lum(end), 'RelTol', 0.05);      % arms matched in LINEAR light
            tc.verifyGreaterThan(lum(128), 3 * max(lum(1), lum(end)));   % neutral centre
        end

        function theScalogramIsCachedPerPageNotPerFrame(tc)
            % It is a whole-page quantity; recomputing it on every frame would dominate the
            % redraw for no benefit.
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'alpha', 'PageIndex', 1);
            E = tc.App.ScalogramData;
            tc.verifySize(E, [tc.App.Page.NumScales, tc.App.Page.NumFrames]);
            tc.App.setFrame(tc.App.Frame + 5);
            tc.verifyEqual(tc.App.ScalogramData, E);            % untouched by scrubbing
            tc.App.onPage(2);
            tc.verifySize(tc.App.ScalogramData, [tc.App.Page.NumScales, tc.App.Page.NumFrames]);
        end

        function advancingPagesKeepsTheBand(tc)
            tc.App = rheome.flowbrowser(tc.B, 'Band', 'beta', 'PageIndex', 1);
            tc.App.onPage(2);
            tc.verifyEqual(tc.App.PageIndex, 2);
            tc.verifyEqual(tc.App.BandName, 'beta');
            tc.verifyEqual(tc.App.Page.PageIndex, 2);
        end

    end
end

% Author: Diellor Basha, 2026
