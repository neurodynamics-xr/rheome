classdef tSelectBursts < matlab.unittest.TestCase
% The two-phase query against a brute-force scan of the same members at full rate:
% zero false negatives, on/off within a sample, censoring at the record edge, and the
% cost report's raw fraction well below one.
%
% Author: Diellor Basha, 2026

    methods (Static)
        function E = brute(fx, b, theta, minDur, chans)
            db = fx.db;  bands = db.bands;  fs = fx.fs;  nT = fx.nT;
            V = db.meta.cfg.VoicesPerOctave;  a = db.meta.cfg.Anchor;
            tfb = rheome.timefilterbank(nT, 'SamplingFrequency', fs, 'VoicesPerOctave', V, 'Anchor', a);
            [H, f] = freqz(tfb);  fc = centerFrequencies(tfb);  kd = kinds(tfb);
            mem = find(strcmp(kd, 'band') & fc >= bands.fLo(b) - 1e-9 & fc < bands.fHi(b) - 1e-9);
            E = [];
            for ch = chans
                X = fft(fx.X(:, ch));  env = zeros(nT, 1);
                for m = mem, h = zeros(nT, 1);  h(1:numel(f)) = H(m, :);  env = max(env, abs(ifft(X .* h))); end
                on = env >= theta;  d = diff([0; on; 0]);  st = find(d == 1);  sp = find(d == -1) - 1;
                for e = 1:numel(st)
                    dur = (sp(e) - st(e) + 1) / fs;
                    if dur >= minDur, E = [E; ch, (st(e)-1)/fs, sp(e)/fs, max(env(st(e):sp(e)))]; end %#ok<AGROW>
                end
            end
        end
    end

    methods (Test)

        function everyBruteForceEventIsFoundWithinASample(tc)
            fx = selFixture();  db = fx.db;
            b = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);
            theta = 1.5e-13;  minDur = 0.3;
            [E, cost] = rheome.select.bursts(db, Band=b, Threshold=theta, MinDuration=minDur, Channels=[1 3]);
            B = tSelectBursts.brute(fx, b, theta, minDur, [1 3]);
            tc.verifyGreaterThan(size(B, 1), 3, 'the fixture must hold events');
            for i = 1:size(B, 1)
                j = find(E.channel_id == B(i, 1) & abs(E.t_on - B(i, 2)) <= 1.5/fx.fs & abs(E.t_off - B(i, 3)) <= 1.5/fx.fs);
                tc.verifyNumElements(j, 1, sprintf('brute event ch %d at %.3f s', B(i, 1), B(i, 2)));
                if isscalar(j), tc.verifyEqual(E.peak(j), B(i, 4), 'RelTol', 1e-3); end   % span vs whole-record evaluation: 5e-5 measured
            end
            tc.verifyEqual(height(E), size(B, 1));                              % and no extra events
            tc.verifyLessThanOrEqual(cost.rawFraction, 1);                       % spans merged: nothing read twice (the fixture is all bursts)
            tc.verifyGreaterThan(cost.rawSamplesRead, 0);
            tc.verifyEqual(cost.phase1.executor, 'frames');
        end

        function aHighThresholdReadsNoRawSamples(tc)
            fx = selFixture();  db = fx.db;
            b = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);
            [E, cost] = rheome.select.bursts(db, Band=b, Threshold=1, MinDuration=0.3);
            tc.verifyEqual(height(E), 0);
            tc.verifyEqual(cost.rawSamplesRead, 0);
        end

        function eventsNearTheRecordEdgeAreCensored(tc)
            fx = selFixture();  db = fx.db;
            b = find(db.bands.fLo <= 3 & 3 < db.bands.fHi);                    % channel 2: 3 Hz throughout
            [E, ~] = rheome.select.bursts(db, Band=b, Threshold=1e-13, MinDuration=0, Channels=2);
            tc.verifyGreaterThan(height(E), 0);
            tc.verifyTrue(any(E.censored));
            tc.verifyTrue(all(E.t_on >= 0 & E.t_off <= fx.nT / fx.fs + 1e-9));
        end

    end
end
