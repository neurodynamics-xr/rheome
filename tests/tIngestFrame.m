classdef tIngestFrame < matlab.unittest.TestCase
% The ingest on the frame bank (the default): integer-octave bands plus the two edge
% bands partition sum x^2 exactly, the sub-band reduce is Parseval-exact over the record,
% roll-up and build are unchanged, and the numbers match the full-rate Morse path where
% the two banks measure the same thing.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 4000
        C  = 2
        X
        cfg
        fb
        bands
        frame
        g
    end

    methods (TestClassSetup)
        function make(tc)
            rng(21);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13*randn(tc.nT, tc.C) + 5e-13*sin(2*pi*10*t) .* (1 + 0.5*sin(2*pi*0.2*t));
            tc.cfg = rheome.ingest.config();                               % Bank = "frame"
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
        end
    end

    methods (Test)

        function theDefaultBankIsTheFrameAndItIsTight(tc)
            tc.verifyClass(tc.fb, 'rheome.timefilterbank');
            tc.verifyEqual(tc.frame.A, 1, 'AbsTol', 1e-12);
            tc.verifyEqual(tc.frame.B, 1, 'AbsTol', 1e-12);
            tc.verifyEqual(tc.frame.bank, 'frame');
        end

        function bandsAreIntegerOctavesOfTheAnchorPlusTwoEdges(tc)
            b = tc.bands;
            tc.verifyEqual(b.kind{1}, 'above');  tc.verifyEqual(b.kind{end}, 'below');
            oc = strcmp(b.kind, 'octave') & ~b.partial;                                 % full octaves
            tc.verifyGreaterThan(nnz(oc), 3);
            tc.verifyEqual(b.fLo(oc), 2.^round(log2(b.fLo(oc))), 'RelTol', 1e-12);   % anchor 1 Hz
            tc.verifyEqual(b.fExtent(oc), ones(nnz(oc), 1), 'AbsTol', 1e-12);
            oc = strcmp(b.kind, 'octave');
            tc.verifyEqual(b.fLo(1:end-1), b.fHi(2:end), 'RelTol', 1e-12);          % contiguous, top first
            all_ = sort(cell2mat(b.scales'));
            tc.verifyEqual(all_, 1:tc.fb.NumMembers);                                 % every member in one band
            tc.verifyTrue(all(cellfun(@numel, b.scales(oc)) <= tc.cfg.VoicesPerOctave));
        end

        function bandEnergiesPartitionTheRecordExactly(tc)
            % sum over bands and tiles of energy == sum x^2 minus the DC and Nyquist bins,
            % per channel, to 1e-10: the frame is tight and nothing is truncated.
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            for c = 1:tc.C
                Xf = fft(tc.X(:, c));
                excluded = (abs(Xf(1))^2 + abs(Xf(tc.nT/2+1))^2) / tc.nT;
                tc.verifyEqual(sum(T0.energy(:, c, :), 'all'), sum(tc.X(:, c).^2) - excluded, 'RelTol', 1e-10);
                tc.verifyEqual(sum(T0.sumX2(:, c)), sum(tc.X(:, c).^2), 'RelTol', 1e-12);
            end
        end

        function perTilePartitionHoldsAtTheNaturalLevelAndAbove(tc)
            % Per tile the identity is a coarse-level statement (design §3.13): at the top
            % level it is exact; at 4 s tiles it holds to a few percent for this signal.
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            T = rheome.ingest.rollup(T0, tc.g);
            top = T{end};
            Xf = fft(tc.X(:, 1));  excluded = (abs(Xf(1))^2 + abs(Xf(tc.nT/2+1))^2) / tc.nT;
            tc.verifyEqual(sum(top.energy(1, 1, :)), top.sumX2(1, 1) - excluded, 'RelTol', 1e-10);
            L4 = T{5};  r = sum(L4.energy(:, 1, :), 3) ./ L4.sumX2(:, 1);
            tc.verifyGreaterThan(median(r), 0.95);  tc.verifyLessThan(median(r), 1.05);
        end

        function aToneSitsInItsOctave(tc)
            t = (0:tc.nT-1)' / tc.fs;  x = sin(2*pi*10*t);
            T0 = rheome.ingest.reduce(x, tc.fb, tc.bands, tc.g, tc.frame);
            E = squeeze(sum(T0.energy, [1 2]));
            b = find(tc.bands.fLo <= 10 & 10 < tc.bands.fHi);
            tc.verifyGreaterThan(E(b) / sum(E), 0.9);                  % one voice (a quarter octave) of overhang into neighbours
        end

        function maskAndCountsBehaveAsOnTheMorsePath(tc)
            mask = true(tc.nT, 1);  mask(1001:1050) = false;           % all of frame 21
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame, Mask=mask);
            tc.verifyEqual(T0.n(21), uint32(0));
            tc.verifyEqual(squeeze(T0.energy(21, :, :)), zeros(tc.C, height(tc.bands)));
            tc.verifyEqual(T0.sumX(21, :), zeros(1, tc.C));
            tc.verifyEqual(sum(double(T0.nCoi(:, end))), tc.nT - 50);  % the "below" band: every VALID sample is inside the cone
        end

        function envelopeMaxBoundsEveryMemberOfTheBand(tc)
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            Cc = wt(tc.fb, tc.X(:, 1));
            for b = 1:height(tc.bands)
                for m = tc.bands.scales{b}
                    fi = min(floor(Cc(m).t / tc.cfg.FrameFloor) + 1, tc.g.K0);
                    per = accumarray(fi, abs(Cc(m).coef), [tc.g.K0 1], @max);
                    tc.verifyTrue(all(single(per) <= T0.envMax(:, 1, b) * (1 + 1e-6)));
                end
            end
        end

        function rollupIdentityHoldsOnTheFramePath(tc)
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            T  = rheome.ingest.rollup(T0, tc.g);
            cfg2 = rheome.ingest.config(FrameFloor = 1);                        % level 2 directly
            g2 = rheome.ingest.grid(tc.nT, tc.fs, cfg2);
            D  = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, g2, tc.frame);
            tc.verifyEqual(T{3}.energy, D.energy, 'RelTol', 1e-12);
            tc.verifyEqual(T{3}.envMax, D.envMax);
            tc.verifyEqual(T{3}.n, D.n);
        end

        function buildWritesAFrameStoreWithoutPaging(tc)
            d = tempname;  mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
            rec = fullfile(d, 'recording.mat');
            F = tc.X.';  Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;              %#ok<NASGU>
            ChannelFlag = [1; 1];  ChannelName = {'A','B'};  ChannelType = {'MEG','MEG'}; %#ok<NASGU>
            nCh = 2;  nT = tc.nT;  Comment = '';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(rec, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
                 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', '-v7.3', '-nocompression');
            f = rheome.ingest.build(rec, rheome.ingest.config(Paged="always"), Verbose=false);   % Paged is ignored for the frame bank
            tc.verifyTrue(endsWith(f, 'ingest__F0.25__V4__frame.mat'));
            m = matfile(f);  meta = m.meta;  bands = m.bands;
            tc.verifyFalse(meta.paged);
            tc.verifyEqual(bands.kind{1}, 'above');
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            gS = m.grid;
            tc.verifyEqual(m.L00_sumX2, T0.sumX2);
            L = find(cellfun(@numel, gS.bandsAt) >= 1, 1) - 1;
            T = rheome.ingest.rollup(T0, tc.g);
            tc.verifyEqual(m.(sprintf('L%02d_energy', L)), T{L+1}.energy(:, :, gS.bandsAt{L+1}));
        end

    end
end
