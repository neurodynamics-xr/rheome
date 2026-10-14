classdef tShowTrack < matlab.unittest.TestCase
% rheome.show.track: the Editor's track course (engine surfaceKit trackValue), value by value.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function profiles(t)
            x = [0 0.25 0.5 1];
            t.verifyEqual(rheome.show.track(struct('Start', 1, 'To', 3, 'Profile', 'ramp'), x), [1 1.5 2 3], 'AbsTol', 1e-12);
            t.verifyEqual(rheome.show.track(struct('Start', 0, 'To', 1, 'Profile', 'ease'), x), [0 (1-cos(pi/4))/2 0.5 1], 'AbsTol', 1e-12);
            t.verifyEqual(rheome.show.track(struct('Start', 0, 'To', 1, 'Profile', 'oscillate'), x), [0 0.5 1 0], 'AbsTol', 1e-12);
            t.verifyEqual(rheome.show.track(struct('Start', 0, 'To', 3, 'Profile', 'step'), x), [0 1 2 3], 'AbsTol', 1e-12);
            t.verifyEqual(rheome.show.track(struct('Start', 2), x), [2 2 2 2]);
        end
        function keys(t)
            K = struct('x', {0, 0.5, 1}, 'v', {0, 1, 0}, 'seg', {'hold', 'linear', ''});
            t.verifyEqual(rheome.show.track(struct('Keys', K), [0.25 0.5 0.75 2]), [0 1 0.5 0], 'AbsTol', 1e-12);
            t.verifyEqual(rheome.show.track(struct('Keys', [0.2 1; 0.8 3]), [0 0.5 1]), [1 2 3], 'AbsTol', 1e-12);
            C = struct('x', {0, 0.3, 0.6, 1}, 'v', {0, 2, 2.2, 5}, 'seg', 'curve');
            v = rheome.show.track(struct('Keys', C), linspace(0, 1, 501));
            t.verifyGreaterThanOrEqual(diff(v), -1e-12);                      % monotone keys: a monotone curve, no overshoot
        end
        function oscillator(t)
            s = struct('Start', 0, 'Osc', struct('rate', 2, 'depth', 0.5, 'phase', 0, 'shape', 'sine'));
            t.verifyEqual(rheome.show.track(s, [0 0.125 0.375]), [0 0.5 -0.5], 'AbsTol', 1e-12);
            s.Osc.shape = 'square';
            t.verifyEqual(rheome.show.track(s, [0.1 0.3]), [0.5 -0.5], 'AbsTol', 1e-12);
        end
        function refuses(t)
            t.verifyError(@() rheome.show.track(struct('Nope', 1), 0), 'rheome:show:track');
            t.verifyError(@() rheome.show.track(struct('To', 1, 'Profile', 'zigzag'), 0), 'rheome:show:track');
        end
    end
end
