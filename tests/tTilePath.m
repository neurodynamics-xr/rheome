classdef tTilePath < matlab.unittest.TestCase
% Track-before-detect on the tile graph (rheome.detect.tilepath) and the picker's frame floor.
%
% What is asserted: the Viterbi path is the exact best local walk (checked against brute force on a
% small graph); it follows a planted walk through frames where the walk is NOT the frame's maximum,
% which is the whole point; it may start and end inside the window; a baseline above the field
% returns no path; extra paths are disjoint; and a still feature is never rolled up below MinFrames.
%
% Author: Diellor Basha, 2026

    properties
        G; nTile
    end

    methods (TestClassSetup)
        function tiles(tc)
            [V, F] = rheome.geom.icosphere(3);
            S = struct('Vertices', V * 0.07, 'Faces', F);
            [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces);
            T = rheome.geom.tree(S, L=L, M=M, MaxDepth=4);
            tc.G = rheome.geom.tiles(T, S, 4, Ruler=rheome.geom.edgegraph(S));
            tc.nTile = numel(tc.G.node_id);
        end
    end

    methods (Test)

        function matchesBruteForceOnASmallWindow(tc)
            rng(2);  nF = 4;  Y = randn(tc.nTile, nF);  b = 0.3;
            out = rheome.detect.tilepath(Y, tc.G, Baseline=b);
            best = tc.brute(Y - b, nF);
            tc.verifyEqual(out.tracks(1).score, best, 'AbsTol', 1e-12);
            x = out.tracks(1);
            tc.verifyEqual(sum(Y(sub2ind(size(Y), x.tiles, x.frames)) - b), best, 'AbsTol', 1e-12);
            st = [x.tiles(1:end-1) x.tiles(2:end)];  st = st(st(:,1) ~= st(:,2), :);
            if ~isempty(st)
                tc.verifyTrue(all(full(tc.G.A(sub2ind(size(tc.G.A), st(:,1), st(:,2)))) > 0));
            end
        end

        function followsAWalkThroughFramesItDoesNotWin(tc)
            % a walk of value 1 along adjacent tiles, with a decoy of 1.3 elsewhere on every third frame:
            % a per-frame argmax jumps to the decoy, the path stays on the walk
            nF = 30;  walk = tc.walk(nF);  Y = zeros(tc.nTile, nF);
            Y(sub2ind(size(Y), walk, (1:nF)')) = 1;
            far = find(~tc.G.A(walk(1), :) & (1:tc.nTile) ~= walk(1));
            far = far(~ismember(far, walk));
            for t = 3:3:nF, Y(far(1), t) = 1.3; end
            out = rheome.detect.tilepath(Y, tc.G, Baseline=0.5);
            tc.verifyEqual(out.tracks(1).tiles, walk);
            [~, am] = max(Y, [], 1);
            tc.verifyGreaterThan(sum(am(:) ~= walk), 5, 'the decoy should win the frame argmax');
        end

        function startsAndEndsInsideTheWindow(tc)
            nF = 40;  walk = tc.walk(10);  Y = -0.2 * ones(tc.nTile, nF);
            Y(sub2ind(size(Y), walk, (16:25)')) = 1;
            x = rheome.detect.tilepath(Y, tc.G).tracks(1);
            tc.verifyEqual([x.firstFrame x.lastFrame], [16 25]);
        end

        function fullSpanRunsFirstToLastFrame(tc)
            nF = 40;  walk = tc.walk(10);  Y = -0.2 * ones(tc.nTile, nF);
            Y(sub2ind(size(Y), walk, (16:25)')) = 1;
            x = rheome.detect.tilepath(Y, tc.G, Span="full").tracks(1);
            tc.verifyEqual([x.firstFrame x.lastFrame], [1 nF]);
            tc.verifyEqual(x.tiles(16:25), walk);                  % and still follows the walk
            st = [x.tiles(1:end-1) x.tiles(2:end)];  st = st(st(:,1) ~= st(:,2), :);
            if ~isempty(st), tc.verifyTrue(all(full(tc.G.A(sub2ind(size(tc.G.A), st(:,1), st(:,2)))) > 0)); end
        end

        function risefallRecoversPlantedTimes(tc)
            fs = 100;  t = (0:599) / fs;  S = @(x) 1 ./ (1 + exp(-x));
            a = 1 + 4 * S((t - 2.0) / 0.12) .* S(-(t - 4.1) / 0.2) + 0.05 * randn(size(t));
            R = rheome.detect.risefall(a, fs, 1.8, 4.3);
            tc.verifyTrue(R.ok);
            tc.verifyEqual([R.tOn R.tOff], [2.0 4.1], 'AbsTol', 0.03);
            tc.verifyEqual([R.tauR R.tauF], [0.12 0.2], 'AbsTol', 0.04);
            tc.verifyGreaterThan(R.r2, 0.95);
            tc.verifyEqual(R.plateau, [R.tOn + 2*R.tauR, R.tOff - 2*R.tauF], 'AbsTol', 1e-12);
            R2 = rheome.detect.risefall(3 * a, fs, 1.8, 4.3);            % height does not move the timing
            tc.verifyEqual([R2.tOn R2.tOff], [R.tOn R.tOff], 'AbsTol', 0.01);
            R3 = rheome.detect.risefall(1e-11 * a, fs, 1.8, 4.3);        % nor does a source-envelope scale
            tc.verifyEqual([R3.tOn R3.tOff R3.tauR], [R.tOn R.tOff R.tauR], 'AbsTol', 0.01);
            tc.verifyEqual(R3.h, 1e-11 * R.h, 'RelTol', 1e-3);
        end

        function aHighBaselineReturnsNoPath(tc)
            Y = rand(tc.nTile, 20);
            out = rheome.detect.tilepath(Y, tc.G, Baseline=2);
            tc.verifyEqual(out.nTracks, 0);
        end

        function extraPathsAreDisjoint(tc)
            rng(4);  Y = randn(tc.nTile, 25);
            out = rheome.detect.tilepath(Y, tc.G, NumPaths=3);
            k = arrayfun(@(x) sub2ind(size(Y), x.tiles, x.frames), out.tracks, 'uni', 0);
            k = vertcat(k{:});
            tc.verifyEqual(numel(unique(k)), numel(k));
            tc.verifyTrue(issorted([out.tracks.score], 'descend'));
        end

        function trackstatsReadsAPath(tc)
            nF = 12;  walk = tc.walk(nF);  Y = zeros(tc.nTile, nF);
            Y(sub2ind(size(Y), walk, (1:nF)')) = 1;
            st = rheome.detect.trackstats(rheome.detect.tilepath(Y, tc.G, Baseline=0.5), tc.G, 10);
            tc.verifyEqual(st.netMM, 1e3 * tc.G.D(walk(1), walk(end)), 'AbsTol', 1e-9);
        end

        function aStillFeatureKeepsMinFrames(tc)
            [V, F] = rheome.geom.icosphere(4);
            S = struct('Vertices', V * 0.07, 'Faces', F);
            [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces);
            Lad = rheome.geom.jointladder(rheome.geom.tree(S, L=L, M=M, MaxDepth=8), Fs=600);
            row = rheome.geom.jointcell(Lad, SigmaMM=30, Hz=10, SpeedMS=0);
            tc.verifyGreaterThanOrEqual(row.framesPerWindow, 8);
            row3 = rheome.geom.jointcell(Lad, SigmaMM=30, Hz=10, SpeedMS=0, MinFrames=3);
            tc.verifyLessThanOrEqual(row3.frameHz, row.frameHz);
        end
    end

    methods
        % A random walk on the tile graph that moves every other frame.
        function w = walk(tc, nF)
            rng(1);  w = zeros(nF, 1);  w(1) = 1;
            for t = 2:nF
                if mod(t, 2), nb = find(tc.G.A(w(t-1), :)); w(t) = nb(randi(numel(nb)));
                else, w(t) = w(t-1); end
            end
        end
        % Best local-walk score by exhaustive enumeration over every start, end and walk.
        function best = brute(tc, Z, nF)
            A = tc.G.A > 0;  best = -Inf;
            for t0 = 1:nF
                for k0 = 1:tc.nTile
                    best = max(best, tc.extend(Z, A, k0, t0, Z(k0, t0), nF));
                end
            end
        end
        function b = extend(tc, Z, A, k, t, s, nF)
            b = s;
            if t == nF, return; end
            for j = [k find(A(k, :))]
                b = max(b, tc.extend(Z, A, j, t + 1, s + Z(j, t + 1), nF));
            end
        end
    end
end

% Author: Diellor Basha, 2026
