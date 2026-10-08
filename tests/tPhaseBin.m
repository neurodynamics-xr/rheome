classdef tPhaseBin < matlab.unittest.TestCase
% rheome.flow.phasebin -- collapse a multirate coefficient series onto a COMMON phase axis.
%
% This is what makes multirate rectangular again. Each filter runs at its own sample rate,
% so the TIME axes disagree; binning by the reference rhythm's phase gives every filter the
% same number of bins per cycle, and the result stacks.
%
% The two accumulators answer DIFFERENT questions and the distinction is the point:
%   complex mean -- keeps only what is PHASE-LOCKED to the reference; the cycle-averaged
%                   flow map, i.e. the thing you animate. Non-locked activity cancels.
%   energy  mean -- keeps ALL power at that phase, locked or not. Cannot cancel.
% Their ratio is a locking measure. A test that only checked one would miss a reader that
% silently computed the other.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function assignsSamplesToBinsByPhase(tc)
            % phase -pi..pi over nBins bins: bin edges are deterministic, so a sample at a
            % known phase must land in a known bin.
            nBins = 4;
            phase = [-pi+0.01, -pi/2+0.01, 0.01, pi/2+0.01];
            A = ones(1, 4);
            out = rheome.flow.phasebin(A, phase, nBins);
            tc.verifyEqual(out.count, [1 1 1 1]);
            tc.verifyEqual(out.phase, (-pi + (0.5:nBins)*(2*pi/nBins)), 'AbsTol', 1e-12);
        end

        function complexMeanKeepsPhaseLockedContentAndCancelsTheRest(tc)
            % Two modes: one rigidly locked to the reference, one with random phase.
            % The locked one must survive the cycle average; the other must not.
            rng(7);
            nT = 6000;  nBins = 30;
            phase = mod(2*pi*(0:nT-1)/40, 2*pi) - pi;         % 40 samples per cycle
            locked   = 2 * exp(1i * phase);                    % rigidly follows the reference
            unlocked = exp(1i * 2*pi*rand(1, nT));             % independent of it
            out = rheome.flow.phasebin([locked; unlocked], phase, nBins);

            tc.verifyGreaterThan(mean(abs(out.mean(1,:))), 1.9);    % locked: survives
            tc.verifyLessThan(   mean(abs(out.mean(2,:))), 0.4);    % unlocked: cancels
        end

        function energyMeanSurvivesWhatTheComplexMeanCancels(tc)
            % The same unlocked mode carries real power. Energy must report it.
            rng(7);
            nT = 6000;  nBins = 30;
            phase = mod(2*pi*(0:nT-1)/40, 2*pi) - pi;
            unlocked = exp(1i * 2*pi*rand(1, nT));
            out = rheome.flow.phasebin(unlocked, phase, nBins);
            tc.verifyLessThan(mean(abs(out.mean)), 0.4);            % cancels
            tc.verifyEqual(mean(out.energy), 1, 'AbsTol', 0.05);    % but the power is there
        end

        function everySampleIsCountedExactlyOnce(tc)
            rng(3);
            nT = 1000;  nBins = 17;
            phase = 2*pi*rand(1,nT) - pi;
            out = rheome.flow.phasebin(randn(5,nT) + 1i*randn(5,nT), phase, nBins);
            tc.verifyEqual(sum(out.count), nT);
        end

        function phaseAtPlusPiWrapsIntoTheFirstBinNotAnExtraOne(tc)
            % The wrap point is exactly where an off-by-one creates a phantom bin.
            out = rheome.flow.phasebin(ones(1,3), [-pi, pi, pi-1e-12], 8);
            tc.verifyLength(out.count, 8);
            tc.verifyEqual(sum(out.count), 3);
            tc.verifyEqual(out.count(1), 2);      % -pi and +pi are the SAME phase
        end

        function emptyBinsAreNaNNotZero(tc)
            % A bin no sample reached is UNKNOWN, not "zero flow". Zero would be averaged
            % into downstream statistics as a real measurement.
            out = rheome.flow.phasebin(ones(1,2), [-pi+0.01, -pi+0.02], 8);
            tc.verifyEqual(out.count(1), 2);
            tc.verifyTrue(all(isnan(out.mean(2:end))));
            tc.verifyTrue(all(isnan(out.energy(2:end))));
        end

        function shapesFollowTheInputModeCount(tc)
            out = rheome.flow.phasebin(randn(13, 500) + 1i*randn(13,500), 2*pi*rand(1,500)-pi, 12);
            tc.verifySize(out.mean,   [13 12]);
            tc.verifySize(out.energy, [13 12]);
            tc.verifySize(out.count,  [1 12]);
        end

        function rejectsMismatchedPhaseLength(tc)
            tc.verifyError(@() rheome.flow.phasebin(ones(3,10), ones(1,9), 8), 'flow:phasebin:size');
        end

    end
end

% Author: Diellor Basha, 2026
