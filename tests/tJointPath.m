classdef tJointPath < matlab.unittest.TestCase
% Tracking across frequency (rheome.detect.jointpath / jointpathnull), no instrument.
%
% What is asserted: the path is the exact best local walk on the (tile x band) product graph (brute force on a
% small window) and every step changes the tile OR the band, never both; a planted chirp ridge -- one place,
% its band climbing -- is followed across bands through frames where a decoy wins, and its sweep is read
% back from the frequency centroid; Norm makes a 1/f field commensurate (no path sinks to the lowest band on
% noise); on noise alone the held-out false rate is near Alpha (two-sided); on one broadband noise through
% overlapping bands (bands share envelope variance) a null with one shift per place shared by all bands holds
% the noise's own pass and cross-band rates; a jump of several bands in one frame is a transit or two paths.
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
            T = rheome.geom.tree(S, L=L, M=M, MaxDepth=3);
            tc.G = rheome.geom.tiles(T, S, 3, Ruler=rheome.geom.edgegraph(S));
            tc.nTile = numel(tc.G.node_id);
        end
    end

    methods (Test)
        function exactOnTheProductGraph(tc)
            rng(4);  nB = 3;  nF = 3;  Y = randn(tc.nTile, nB, nF);  b = 0.4;
            out = rheome.detect.jointpath(Y, tc.G, Baseline=b);
            A = tc.G.A > 0;  A(1:tc.nTile+1:end) = false;
            n = tc.nTile*nB;  Z = reshape(Y, n, nF) - b;
            nb = @(j) tc.joint(j, A, nB);
            best = -Inf;                                               % brute force: every walk of 1..nF frames
            for t0 = 1:nF
                for j = 1:n, best = max(best, tc.walk(Z, j, t0, nF, nb)); end
            end
            tc.verifyEqual(out.tracks(1).score, best, 'AbsTol', 1e-12);
            x = out.tracks(1);
            dk = diff(x.tiles) ~= 0;  db = diff(x.bands) ~= 0;
            tc.verifyFalse(any(dk & db), 'a step changes the tile or the band, not both');
            tc.verifyTrue(all(abs(diff(x.bands)) <= 1));
        end

        function followsAChirpAcrossBands(tc)
            % one tile, the band climbing 1 -> 7 over 42 frames (6 frames per band); a decoy elsewhere wins every
            % 4th frame; the walk must stay on the ridge and read a sweep of +0.5 oct per 6 frames
            nB = 7;  nF = 42;  rate = 12;  lf = 2:0.5:5;  k0 = 5;  Y = zeros(tc.nTile, nB, nF);
            bt = min(floor((0:nF-1)/6) + 1, nB);
            for t = 1:nF
                Y(k0, bt(t), t) = 1;
                if bt(t) > 1, Y(k0, bt(t) - 1, t) = 0.3; end            % a band's neighbours see some of it
                if bt(t) < nB, Y(k0, bt(t) + 1, t) = 0.3; end
            end
            far = find(~tc.G.A(k0, :) & (1:tc.nTile) ~= k0, 1);
            for t = 4:4:nF, Y(far, 4, t) = 1.3; end
            out = rheome.detect.jointpath(Y, tc.G, Baseline=0.5, LogF=lf, SampleRate=rate);
            x = out.tracks(1);
            tc.verifyEqual(unique(x.tiles), k0);
            tc.verifyEqual(x.bands(:), bt(x.frames)');
            tc.verifyEqual([x.firstFrame x.lastFrame], [1 nF]);
            tc.verifyEqual(x.sweepOctS, 0.5*rate/6, 'RelTol', 0.15);      % 1 oct/s planted
            tc.verifyGreaterThanOrEqual(x.octSpan, 2.5);
            ts = rheome.detect.trackstats(out, tc.G, rate);               % the tile walk is still a trackstats input
            tc.verifyEqual(ts.netMM, 0);
            tc.verifyGreaterThan(x.netOct, 1.5);                          % first third bands 1-3, last third 5-7
            tc.verifyEqual(x.jumpOct, 0.5, 'AbsTol', 1e-9);
        end

        function aStillToneReadsNoSweep(tc)
            nB = 7;  nF = 40;  lf = 2:0.5:5;  Y = 0.1*rand(tc.nTile, nB, nF);  Y(3, 4, :) = 1;  Y(3, [3 5], :) = 0.3;
            x = rheome.detect.jointpath(Y, tc.G, Baseline=0.5, LogF=lf, SampleRate=10).tracks(1);
            tc.verifyEqual(unique(x.bands), 4);
            tc.verifyEqual(x.sweepOctS, 0, 'AbsTol', 1e-9);
            tc.verifyEqual(x.logf(1), 3.5, 'AbsTol', 1e-9);              % symmetric neighbours: the centre
            tc.verifyEqual([x.netOct x.jumpOct], [0 0], 'AbsTol', 1e-9);
        end

        function nullNormalisesAndHoldsTheFalseRate(tc)
            % 1/f noise: band b is (8/2^b) times louder. Without Norm the best path sits in band 1; with the null's
            % Norm the bands are used about evenly, and the held-out false rate stays near Alpha
            rng(11);  nB = 7;  nF = 20;  rate = 10;  nWin = 400;  lf = 2:0.5:5;
            sc = reshape(8 ./ 2.^(1:nB), 1, nB);
            Ys = sc .* abs(filter(1, [1 -0.6], randn(tc.nTile, nB, nF*nWin), [], 3));
            N = rheome.detect.jointpathnull(Ys, tc.G, nF, rate, LogF=lf, Alpha=0.1);
            tc.verifySize(N.norm, [tc.nTile nB]);
            tc.verifyGreaterThan(N.falseTestS, 0.03);                    % alpha 0.1/s; binomial spread at 200 windows
            tc.verifyLessThan(N.falseTestS, 0.25);
            tc.verifyLessThanOrEqual(N.crossTestS, N.falseTestS);
            tc.verifyTrue(isfinite(N.crossOct) && N.crossOct > 0, 'the cross criterion comes from the null');
            tc.verifyEqual(N.falseTestW, N.falseTestS*N.windowS, 'AbsTol', 1e-12);
            tc.verifyEqual(sum(N.bandOcc), 1, 'AbsTol', 1e-12);
            raw = zeros(1, nB);  nor = raw;
            for w = 1:60
                Yw = Ys(:, :, (w-1)*nF + (1:nF));
                p = rheome.detect.jointpath(Yw, tc.G, Baseline=quantile(Yw(:), 0.9));  raw(p.tracks(1).bands) = raw(p.tracks(1).bands) + 1;
                p = rheome.detect.jointpath(Yw, tc.G, Null=N);  nor(p.tracks(1).bands) = nor(p.tracks(1).bands) + 1;
            end
            tc.verifyGreaterThan(raw(1)/sum(raw), 0.9, 'raw 1/f values sink to the lowest band');
            tc.verifyLessThan(max(nor)/sum(nor), 0.5, 'normalised, no band dominates');
        end

        function sharedEnvelopeBandsHoldTheNull(tc)
            % one broadband noise per tile through octave band-passes at half-octave steps (adjacent envelopes
            % correlate); the null shifts each tile by ONE shift shared by every band (measure_freqtrack's policy,
            % per LB mode there). Every passing path in the noise is false: its pass and cross rates must lie
            % within 3 SD of the null's held-out rates.
            rng(7);  fs = 150;  bs = 8;  rate = fs/bs;  nF = floor(2*rate);  nWin = 400;  lf = 2:0.5:5;  nB = numel(lf);
            nS = nF*bs*nWin;  x = randn(nS, tc.nTile);  Y = zeros(tc.nTile, nB, nF*nWin);
            for b = 1:nB
                [z, p, k] = butter(4, 2.^(lf(b) + [-0.5 0.5])/(fs/2), 'bandpass');  [sos, g] = zp2sos(z, p, k);
                e = abs(hilbert(filtfilt(sos, g, x)));
                Y(:, b, :) = reshape(squeeze(mean(reshape(e, bs, [], tc.nTile), 1))', tc.nTile, 1, []);
            end
            sh = randi([nF, nF*(nWin - 1)], tc.nTile, 1);  Ys = Y;
            for k = 1:tc.nTile, Ys(k, :, :) = circshift(Y(k, :, :), sh(k), 3); end
            N = rheome.detect.jointpathnull(Ys, tc.G, nF, rate, LogF=lf, Alpha=0.1);
            ps = false(nWin, 1);  px = ps;
            for w = 1:nWin
                q = rheome.detect.jointpath(Y(:, :, (w-1)*nF + (1:nF)), tc.G, Null=N, LogF=lf, SampleRate=rate);
                if q.nTracks, ps(w) = q.tracks(1).score > N.threshold;  px(w) = ps(w) && abs(q.tracks(1).netOct) >= N.crossOct; end
            end
            nT = numel(N.testScores);
            for c = {mean(ps), N.falseTestW, 'pass'; mean(px), N.crossTestW, 'cross'}'
                p0 = max(c{2}, 1/nT);  sd = sqrt(p0*(1 - p0)*(1/nWin + 1/nT));
                tc.verifyLessThan(abs(c{1} - c{2}), 3*sd, sprintf('%s rate: noise %.3f, null %.3f', c{3}, c{1}, c{2}));
            end
        end

        function aJumpOfSeveralBands(tc)
            % band 1 then band 5 at one tile, no frame between: long ridges are joined by a transit at one band per
            % frame (3 frames at value 0); short ones are two paths
            nB = 7;  lf = 2:0.5:5;  k0 = 2;  Y = zeros(tc.nTile, nB, 20);  Y(k0, 1, 1:10) = 1;  Y(k0, 5, 11:20) = 1;
            x = rheome.detect.jointpath(Y, tc.G, Baseline=0.5, LogF=lf, SampleRate=10).tracks(1);
            tc.verifyEqual(unique(x.tiles), k0);
            tc.verifyEqual([x.frames(1) x.frames(end) x.bands(1) x.bands(end)], [1 20 1 5]);
            tc.verifyTrue(all(ismember(diff(x.bands), [0 1])), 'the transit climbs one band per frame');
            tc.verifyEqual(sum(x.values == 0), 3);
            tc.verifyLessThanOrEqual(x.jumpOct, 1 + 1e-9, 'the centroid steps <= 2 bands (it sees a ridge next door)');
            tc.verifyEqual(x.netOct, 2, 'AbsTol', 1e-9);
            Y = zeros(tc.nTile, nB, 4);  Y(k0, 1, 1:2) = 1;  Y(k0, 5, 3:4) = 1;
            out = rheome.detect.jointpath(Y, tc.G, Baseline=0.5, LogF=lf, NumPaths=2);
            tc.verifyEqual(out.nTracks, 2);
            tc.verifyEqual(sort([unique(out.tracks(1).bands) unique(out.tracks(2).bands)]), [1 5]);
        end

        function emptyWhenNothingBeatsTheBaseline(tc)
            out = rheome.detect.jointpath(zeros(tc.nTile, 3, 5), tc.G, Baseline=1);
            tc.verifyEqual(out.nTracks, 0);
            tc.verifyTrue(all(isfield(out.tracks, ["bands" "logf" "sweepOctS" "octSpan" "netOct" "jumpOct"])));
        end
    end

    methods
        function J = joint(tc, j, A, nB) %#ok<INUSD>
            k = mod(j - 1, tc.nTile) + 1;  b = (j - k)/tc.nTile + 1;
            J = j;
            J = [J; find(A(:, k)) + tc.nTile*(b - 1)];
            if b > 1, J(end+1) = j - tc.nTile; end
            if b < nB, J(end+1) = j + tc.nTile; end
        end

        function s = walk(tc, Z, j, t, nF, nb)
            % best score of a walk that starts at node j, frame t, and may stop at any frame
            s = Z(j, t);
            if t < nF
                c = -Inf;
                for i = nb(j)', c = max(c, tc.walk(Z, i, t + 1, nF, nb)); end
                s = s + max(c, 0);
            end
        end
    end
end
