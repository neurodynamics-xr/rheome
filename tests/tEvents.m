classdef tEvents < matlab.unittest.TestCase
% Event table (risefall on an event's own area-weighted envelope) and event-locked time-frequency power.

    methods (Test)
        function eventtableRecoversPlantedBoundaries(tc)
            rng(3);  fs = 300;  t = (0:fs*12-1)/fs;  Sg = @(x) 1 ./ (1 + exp(-x));
            ep = Sg((t - 5)/0.15) .* Sg(-(t - 6.4)/0.25);
            Y = 0.05*randn(4, numel(t)) + 1;  Y(2, :) = Y(2, :) + 4*ep;  Y(3, :) = Y(3, :) + 2*ep;
            A = struct('L', [1; 1; 3; 1]);
            Ev = table("L", {[2 3]}, 5.1, 6.2, 'VariableNames', {'hemi','tiles','t0','t1'});
            T = rheome.detect.eventtable(Ev, struct('L', Y), A, fs);
            tc.verifyTrue(T.fit);
            tc.verifyEqual([T.tOn T.tOff], [5 6.4], 'AbsTol', 0.03);
            tc.verifyEqual([T.tauR T.tauF], [0.15 0.25], 'AbsTol', 0.05);
            tc.verifyEqual(T.areaMM2, 4e6);  tc.verifyEqual(T.nTiles, 2);
            tc.verifyEqual(T.height, (4*1 + 2*3)/4, 'RelTol', 0.1);                 % area-weighted mean of the tiles
            Ev.t0 = Ev.t0 + 10;  Ev.t1 = Ev.t1 + 10;                                % Y's first sample at 10 s
            T2 = rheome.detect.eventtable(Ev, struct('L', Y), A, fs, T0=10);
            tc.verifyEqual(T2.tOn, T.tOn + 10, 'AbsTol', 1e-6);
        end

        function noEpisodeIsNotAFit(tc)
            rng(5);  fs = 300;  Y = struct('L', 1 + 0.2*randn(2, fs*8));
            Ev = table("L", {1}, 3, 4, 'VariableNames', {'hemi','tiles','t0','t1'});
            T = rheome.detect.eventtable(Ev, Y, struct('L', [1; 1]), fs);
            tc.verifyFalse(T.fit);
        end

        function tauBoundIsFlagged(tc)
            fs = 300;  t = (0:fs*10-1)/fs;  Sg = @(x) 1 ./ (1 + exp(-x));
            a = 1 + 3*Sg((t - 4)/1.5) .* Sg(-(t - 6)/0.1);
            R = rheome.detect.risefall(a, fs, 4, 6, TauMax=0.5);
            tc.verifyTrue(R.atBound);
            R2 = rheome.detect.risefall(a, fs, 4, 6, TauMax=3);
            tc.verifyFalse(R2.atBound);
        end

        function quadratureReadsAnEpisodeOnABackground(tc)
            % an episode on an incoherent background adds in power: the linear model reads its flanks as sharper
            rng(2);  fs = 300;  t = (0:fs*8-1)/fs;  Sg = @(x) 1 ./ (1 + exp(-x));
            a = sqrt(1 + (3*Sg((t - 3)/0.3) .* Sg(-(t - 5)/0.3)).^2) + 0.02*randn(size(t));
            Rq = rheome.detect.risefall(a, fs, 3, 5, Model="quadrature");
            Rl = rheome.detect.risefall(a, fs, 3, 5);
            tc.verifyEqual([Rq.tauR Rq.tauF], [0.3 0.3], 'AbsTol', 0.04);
            tc.verifyEqual([Rq.tOn Rq.tOff], [3 5], 'AbsTol', 0.03);
            tc.verifyLessThan(Rl.tauR, 0.9*0.3);  tc.verifyGreaterThan(Rl.tOn - 3, 0.05);   % linear: sharper, late
        end

        function posteriorReadsAnEpisodeOnABeatingBackground(tc)
            % the episode beats with an ongoing rhythm of its own band: |p + n|, n a slowly drifting complex
            % background (episode 20 dB over its median: the level at which the through-MEG plants read tau at x0.97-1.02); over planted episodes the posterior medians stay within a
            % time tile and tau within 0.8-1.25 (the through-MEG check at 10 dB is measure_eventplants')
            fs = 300;  t = (0:fs*8-1)/fs;  Sg = @(x) 1 ./ (1 + exp(-x));  e = zeros(20, 4);  r = e;
            for k = 1:20
                rng(k);  on = 3 + 0.5*rand;  off = on + 1 + rand;  tau = 0.1 + 0.3*rand(1, 2);
                n = movmean(randn(1, numel(t)) + 1i*randn(1, numel(t)), round(0.2*fs));  n = n / median(abs(n));
                a = abs(10*Sg((t - on)/tau(1)) .* Sg(-(t - off)/tau(2)) + n);
                R = rheome.detect.risefall(a, fs, on + 0.1, off - 0.1, Model="posterior");
                tc.verifyTrue(R.ok);
                if k <= 5                                            % capped grid: coarse, then 25 ms over the support
                    Rc = rheome.detect.risefall(a, fs, on + 0.1, off - 0.1, Model="posterior", GridMax=60, FramesMax=100);
                    tc.verifyTrue(Rc.ok);  tc.verifyEqual([Rc.tOn Rc.tOff], [R.tOn R.tOff], 'AbsTol', 0.025);
                end
                e(k, :) = [R.tOn - on, R.tOff - off, R.tauR / tau(1), R.tauF / tau(2)];
            end
            tc.verifyLessThan(median(abs(e(:, 1:2))), [0.165 0.165]);
            tc.verifyGreaterThan(median(e(:, 3:4)), [0.8 0.8]);  tc.verifyLessThan(median(e(:, 3:4)), [1.25 1.25]);
        end

        function posteriorTimeIsBoundedOnALongEvent(tc)
            % a 120-s episode: the capped grid and frames keep the fit to seconds, and it still finds the edges
            fs = 300;  t = (0:fs*126-1)/fs;  Sg = @(x) 1 ./ (1 + exp(-x));  rng(7);
            n = movmean(randn(1, numel(t)) + 1i*randn(1, numel(t)), round(0.2*fs));  n = n / median(abs(n));
            a = abs(10*Sg((t - 3)/0.2) .* Sg(-(t - 123)/0.3) + n);
            tic;  R = rheome.detect.risefall(a, fs, 3.1, 122.9, Model="posterior");  s = toc;
            tc.verifyTrue(R.ok);  tc.verifyEqual([R.tOn R.tOff], [3 123], 'AbsTol', 0.165);
            tc.verifyLessThan(s, 120);
        end

        function eventlockedFindsTheBurstFrequency(tc)
            rng(1);  fs = 250;  t = (0:fs*20-1)/fs;  x = 0.1*randn(2, numel(t));
            on = 8;  b = t >= on & t < on + 1;  x(1, b) = x(1, b) + sin(2*pi*11*t(b));
            P = rheome.detect.eventlocked(x, fs, [on; 0.5], [1 0], Freqs=(2:30)');
            tc.verifySize(P.P, [29 numel(P.t) 2]);
            tc.verifyEqual(P.t([1 end]), [-1 1], 'AbsTol', 1e-9);
            tc.verifyTrue(all(isnan(P.P(:, :, 2)), 'all'));                         % segment leaves the data
            m = mean(P.P(:, P.t > 0.2 & P.t < 0.8, 1), 2);  [~, k] = max(m);
            tc.verifyEqual(P.f(k), 11);
            tc.verifyGreaterThan(mean(P.P(10, P.t > 0.2 & P.t < 0.8, 1)), 20*mean(P.P(10, P.t < -0.4, 1)));
        end
    end
end
