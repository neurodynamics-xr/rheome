classdef tSpindle < matlab.unittest.TestCase
% rheome.detect.spindle -- sleep spindle detection on band envelope, with duration gating.
%
% Tested against synthetic bursts of KNOWN onset, duration and frequency, because the
% failure mode on real data is silent: a detector that fires on eye blinks still returns
% plausible-looking events with plausible durations, and only the topography gives it away.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200; nT = 200*60;          % one minute
    end

    methods
        function [x, t0, dur] = burst(tc, tStart, dur, f, amp)
            rng(9);
            x = 0.1*randn(1, tc.nT);
            t = (0:tc.nT-1)/tc.fs;
            i = t >= tStart & t < tStart+dur;
            env = zeros(1, tc.nT);
            env(i) = amp * sin(pi*(t(i)-tStart)/dur).^2;      % smooth on/off
            x = x + env .* sin(2*pi*f*t);
            t0 = tStart;
        end
    end

    methods (Test)

        function findsASyntheticSpindleAtTheRightTime(tc)
            [x, t0, dur] = tc.burst(20, 1.0, 13, 3);
            E = rheome.detect.spindle(x, tc.fs);
            tc.verifyEqual(height(E), 1);
            tc.verifyEqual(E.startSec, t0, 'AbsTol', 0.25);
            tc.verifyEqual(E.durationSec, dur, 'AbsTol', 0.4);
            tc.verifyEqual(E.peakFreq, 13, 'AbsTol', 1.5);
        end

        function rejectsBurstsShorterThanTheMinimumDuration(tc)
            % A 0.2 s burst is not a spindle; without duration gating every envelope
            % excursion qualifies and the count becomes a threshold artefact.
            %
            % ⚠ AMPLITUDE MUST BE REALISTIC HERE, and that is a property of the detector
            % rather than of the fixture. Smoothing the envelope and extending each event
            % down to the boundary threshold WIDENS it, and the wider the louder: a 0.2 s
            % burst at 109 robust-SD above baseline stays above the boundary for 2.4 s and
            % is reported as a 0.80 s event. Envelope detectors therefore over-report the
            % duration of very high-amplitude short events, and a duration read off one
            % should not be trusted near the gate.
            x = tc.burst(20, 0.2, 13, 0.35);          % ~2.5x background, as a spindle is
            tc.verifyEmpty(rheome.detect.spindle(x, tc.fs), 'a 0.2 s burst must be rejected');
        end

        function rejectsBurstsLongerThanTheMaximumDuration(tc)
            % Sustained alpha or an artefact can hold the envelope up for many seconds.
            x = tc.burst(10, 8, 13, 3);
            E = rheome.detect.spindle(x, tc.fs);
            tc.verifyTrue(isempty(E) || all(E.durationSec <= 3.01));
        end

        function findsSeveralSeparateSpindles(tc)
            x = tc.burst(10, 1, 13, 3) + tc.burst(30, 1.2, 12, 3) + tc.burst(45, 0.8, 14, 3);
            E = rheome.detect.spindle(x, tc.fs);
            tc.verifyGreaterThanOrEqual(height(E), 3);
            tc.verifyEqual(sort(round(E.startSec)).', [10 30 45], 'AbsTol', 1);
        end

        function theFalsePositiveRateOnNoiseIsZeroAtTheDefaults(tc)
            % ⚠ NOT free. The envelope of narrowband noise is Rayleigh-distributed, so it is
            % heavy-tailed relative to its MAD and a naive threshold fires steadily on
            % background. What controls this is the BOUNDARY threshold, because it sets how
            % long a noise excursion is extended and therefore whether it clears the
            % duration gate -- measured, loSD 1.5 gives 0.4 false events/min where loSD 2.5
            % gives none, at identical sensitivity.
            rng(1);
            n = 0.1*randn(1, tc.nT);                  % the SAME noise for both settings,
            tc.verifyEmpty(rheome.detect.spindle(n, tc.fs)); % or this compares draws, not settings
            loose = rheome.detect.spindle(n, tc.fs, struct('loSD',1.0));
            tc.verifyGreaterThan(height(loose), 0, 'a loose boundary must let noise through');
        end

        function theThresholdBaselineIsGLOBALNotPerWindow(tc)
            % ⚠ THE DESIGN DECISION THAT DECIDES WHETHER STAGE CONTRAST SURVIVES.
            % Thresholding each segment against ITSELF finds "spindles" in wake too, because
            % wake's own envelope percentile is low. The threshold must come from a shared
            % baseline, so a quiet segment yields nothing.
            loud  = tc.burst(20, 1, 13, 3);
            quiet = 0.1*randn(1, tc.nT);
            both  = [loud quiet];
            E = rheome.detect.spindle(both, tc.fs);
            inQuiet = E.startSec > tc.nT/tc.fs;
            tc.verifyFalse(any(inQuiet), 'the quiet half must yield no detections');
        end

        function multichannelReportsThePerChannelEvents(tc)
            x1 = tc.burst(20, 1, 13, 3);
            x2 = 0.1*randn(1, tc.nT);
            E = rheome.detect.spindle([x1; x2], tc.fs);
            tc.verifyTrue(all(E.channel == 1));
            tc.verifyGreaterThanOrEqual(height(E), 1);
        end

        function timesAreReportedInTheORIGINALSampleIndex(tc)
            % The detector decimates internally for speed and conditioning; if it reported
            % the decimated index every downstream alignment would be silently wrong.
            [x, t0] = tc.burst(20, 1, 13, 3);
            E = rheome.detect.spindle(x, tc.fs);
            tc.verifyEqual(E.startSample(1)/tc.fs, t0, 'AbsTol', 0.25);
            tc.verifyLessThanOrEqual(E.endSample(1), numel(x));
        end

        function rejectsAnImpossibleBand(tc)
            tc.verifyError(@() rheome.detect.spindle(randn(1,1000), 20), 'detect:spindle:band');
        end

    end
end

% Author: Diellor Basha, 2026
