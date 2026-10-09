classdef tEventSensors < matlab.unittest.TestCase
% A tracked event labels the sensor data: the samples of its frames and the channels its tile drives.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function theEventsChannelsAreTheOnesItsTileDrives(tc)
            [track, G, Gain, K, X, fs] = i_case();
            E = rheome.detect.eventsensors(track, G, Gain, K, X, FrameSamples=fs, Fraction=0.5);
            tc.verifyEqual(E.tiles, [2; 3]);
            tc.verifyEqual(sort(E.channels{1}), [1; 2]);            % tile 2 drives channels 1-2
            tc.verifyEqual(sort(E.channels{2}), [3; 4]);            % tile 3 drives channels 3-4
            tc.verifyTrue(all(E.mask([1 2], 11:20), 'all'));
            tc.verifyFalse(any(E.mask(:, [1:10 31:end]), 'all'));   % only the event's frames
            tc.verifyFalse(any(E.mask([3 4], 11:20), 'all'));
        end

        function theEventSignalIsTheTilesEstimateForwarded(tc)
            [track, G, Gain, K, X, fs] = i_case();
            E = rheome.detect.eventsensors(track, G, Gain, K, X, FrameSamples=fs);
            r = reshape((find(G.P(:, 2))' - 1) * 3 + (1:3)', [], 1);
            tc.verifyEqual(E.Xevent(:, 11:20), Gain(:, r) * K(r, :) * X(:, 11:20), 'AbsTol', 1e-12);
            tc.verifyEqual(E.Xevent(:, 1:10), zeros(6, 10));
        end

        function mismatchedShapesAreAnError(tc)
            [track, G, Gain, K, X, fs] = i_case();
            tc.verifyError(@() rheome.detect.eventsensors(track, G, Gain(1:5, :), K, X, FrameSamples=fs), ...
                           'detect:eventsensors:size');
        end

        function theFigureDrawsTracesAndTopography(tc)
            [track, G, Gain, K, X, fs] = i_case();
            E = rheome.detect.eventsensors(track, G, Gain, K, X, FrameSamples=fs);
            [h, info] = rheome.show.eventsensors(X, (0:size(X,2)-1)/100, E, randn(6, 2), 'NumChannels', 4, 'Visible', 'off');
            tc.addTeardown(@close, h);
            tc.verifyEqual(sort(info.carrying), (1:4)');
            tc.verifyNumElements(findobj(h, 'Type', 'axes'), 2);
        end
    end
end

% 4 tiles x 2 vertices (free orientation, 24 source rows), 6 channels; tile k drives channels 2k-3, 2k-2
% (k = 2, 3) and nothing else. The track sits on tile 2 in frame 2, tile 3 in frame 3; 10 samples/frame.
function [track, G, Gain, K, X, fs] = i_case()
    nV = 8;  tileOf = repelem((1:4)', 2);
    G = struct('P', sparse(1:nV, tileOf, true, nV, 4), 'tileOf', tileOf);
    Gain = 0.01 * ones(6, 3*nV);
    for k = 2:3, r = (find(tileOf == k)' - 1) * 3 + (1:3)';  Gain(2*k-3:2*k-2, r(:)) = 1; end
    K = Gain';  rng(3);  X = randn(6, 40);
    fs = [1 10; 11 20; 21 30; 31 40];
    track = struct('frames', [2 3], 'tiles', [2 3]);
end

% Author: Diellor Basha, 2026
