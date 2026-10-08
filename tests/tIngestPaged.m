classdef tIngestPaged < matlab.unittest.TestCase
% The paged reduce against the whole-record reduce on the same signal. Counts and all-pass
% statistics must be identical; band statistics differ by the wavelet tail outside the
% halo, which is measured here and pinned at the level the notes record.
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
        whole
    end

    methods (TestClassSetup)
        function makeSignal(tc)
            rng(13);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13*randn(tc.nT, tc.C) + 5e-13*sin(2*pi*10*t) .* (1 + 0.5*sin(2*pi*0.2*t));
            tc.cfg = rheome.ingest.config(Bank="morse", MinPageLength=1, PageOverhead=0.5);     % 1 s pages where allowed
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            tc.whole = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
        end
    end

    methods

        function T0 = paged(tc, varargin)
            plan = rheome.ingest.pages(tc.g, tc.bands, tc.cfg);
            readfcn = @(a, b) tc.X(a:b, :);
            T0 = rheome.ingest.reducepaged(readfcn, tc.C, tc.bands, tc.g, tc.frame, plan, varargin{:});
        end

    end

    methods (Test)

        function countsAndAllPassStatisticsAreIdentical(tc)
            P = tc.paged();
            W = tc.whole;
            tc.verifyEqual(P.n, W.n);
            tc.verifyEqual(P.nCoi, W.nCoi);
            tc.verifyEqual(P.sumX, W.sumX);
            tc.verifyEqual(P.sumX2, W.sumX2);
            tc.verifyEqual(P.absMax, W.absMax);
            tc.verifyEqual(P.min, W.min);
            tc.verifyEqual(P.max, W.max);
        end

        function outsideTheConePagedAndWholeAgreeAndInsideTheStoreSaysSo(tc)
            % The strong statement of design §3.12: a paged tile equals the whole-record
            % tile to single precision unless the record edge lies within the band's
            % support of it -- and those tiles carry nCoi > 0. Bands at Lmax run whole
            % in both paths and are identical.
            % ⚠ Compared with RelTol, never with ./max(x, eps): the energies here are
            % ~1e-26, far below eps, and max(x, eps) would report a spurious zero.
            P = tc.paged();  W = tc.whole;
            plan = rheome.ingest.pages(tc.g, tc.bands, tc.cfg);
            for r = 1:height(plan)
                for b = plan.bands{r}
                    if plan.level(r) == tc.g.Lmax
                        tc.verifyEqual(P.energy(:, :, b), W.energy(:, :, b), sprintf('band %d whole', b));
                        tc.verifyEqual(P.envMax(:, :, b), W.envMax(:, :, b), sprintf('band %d whole', b));
                    else
                        clear_ = W.nCoi(:, b) == 0;
                        tc.verifyTrue(any(clear_), sprintf('band %d has tiles outside the cone', b));
                        tc.verifyEqual(P.energy(clear_, :, b), W.energy(clear_, :, b), 'RelTol', 1e-5, sprintf('band %d energy outside the cone', b));
                        tc.verifyEqual(double(P.envMax(clear_, :, b)), double(W.envMax(clear_, :, b)), 'RelTol', 1e-5, sprintf('band %d envMax outside the cone', b));   % double: RelTol is ignored on single
                        inside = ~clear_;
                        if any(inside)
                            tc.verifyTrue(all(W.energy(inside, :, b) > 0, 'all'));   % the cone tiles are real tiles
                        end
                    end
                end
            end
        end

        function aMaskIsHonouredIdentically(tc)
            mask = true(tc.nT, 1);  mask(701:760) = false;
            P = tc.paged(Mask=mask);
            W = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame, Mask=mask);
            tc.verifyEqual(P.n, W.n);
            tc.verifyEqual(P.nCoi, W.nCoi);
            tc.verifyEqual(P.sumX, W.sumX);
            tc.verifyEqual(P.min, W.min);
            tc.verifyEqual(P.energy(15, :, :), W.energy(15, :, :));     % a fully masked frame: zeros
            tc.verifyEqual(P.envMax(15, :, :), W.envMax(15, :, :));
        end

        function pagesTileTheRecordExactlyOnce(tc)
            % Every core sample is read exactly once per job; spans overlap only by halos.
            plan = rheome.ingest.pages(tc.g, tc.bands, tc.cfg);
            for r = 1:height(plan)
                L = plan.level(r);  n = tc.g.F * 2^L;
                cores = zeros(tc.nT, 1);
                for p = 1:plan.nPages(r)
                    a = (p-1)*n + 1;  b = min(p*n, tc.nT);
                    cores(a:b) = cores(a:b) + 1;
                end
                tc.verifyTrue(all(cores == 1));
            end
        end

        function theBuildTakesThePagedPathWhenForcedAndTagsTheFile(tc)
            d = tempname;  mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
            rec = fullfile(d, 'recording.mat');
            F = tc.X.';  Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;              %#ok<NASGU>
            ChannelFlag = [1; 1];  ChannelName = {'A','B'};  ChannelType = {'MEG','MEG'}; %#ok<NASGU>
            nCh = 2;  nT = tc.nT;  Comment = '';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(rec, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
                 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', '-v7.3', '-nocompression');
            cfgP = rheome.ingest.config(Bank="morse", MinPageLength=1, PageOverhead=0.5, Paged="always");
            fP = rheome.ingest.build(rec, cfgP, Verbose=false);
            tc.verifyTrue(endsWith(fP, '__paged.mat'));
            m = matfile(fP);  meta = m.meta;
            tc.verifyTrue(meta.paged);
            tc.verifyGreaterThan(height(meta.plan), 1);
            P = tc.paged();
            gP = m.grid;
            tc.verifyEqual(m.L00_energy, P.energy(:, :, gP.bandsAt{1}));   % whole variable: a single band is 2-D on disk
            tc.verifyEqual(m.L00_n, P.n);
            fW = rheome.ingest.build(rec, tc.cfg, Verbose=false);            % auto -> whole here
            tc.verifyFalse(endsWith(fW, '__paged.mat'));
            mW = matfile(fW);  metaW = mW.meta;
            tc.verifyFalse(metaW.paged);
        end

    end
end
