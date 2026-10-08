classdef tFlowFeatures < matlab.unittest.TestCase
% @flowfeatures -- windowed, band- and scale-stratified tabulation of the flow tensor.
%
% The correctness anchor is that the Gram route equals direct per-frame computation. Every
% feature here is a quadratic form, and a quadratic form commutes with summation over a
% window, so the [Ks x Ks] coefficient covariance is a SUFFICIENT STATISTIC -- but only if
% the algebra is right, and a wrong-but-plausible energy is invisible downstream.
%
% Author: Diellor Basha, 2026

    properties
        B; C = 8; V = 60; Ks = 10; nT = 900; fs = 400; nROI = 4;
    end

    methods (TestClassSetup)
        function makeBundle(tc)
            d = tempname; mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            f = fullfile(d, 'recording.mat');
            rng(21);
            F = randn(tc.C, tc.nT);                                    %#ok<NASGU>
            Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;                  %#ok<NASGU>
            ChannelFlag = ones(tc.C,1);                                %#ok<NASGU>
            ChannelName = compose('C%d', 1:tc.C);                      %#ok<NASGU>
            ChannelType = repmat({'MEG'},1,tc.C);                      %#ok<NASGU>
            nCh = tc.C; nT = tc.nT; Comment=''; source='';             %#ok<NASGU>
            save(f,'F','Time','sfreq','ChannelFlag','ChannelName','ChannelType', ...
                 'nCh','nT','Comment','source','-v7.3','-nocompression');

            [Phi,~] = qr(randn(tc.V, tc.Ks), 0);
            Lambda  = linspace(1,400,tc.Ks)';
            wV      = 0.5 + rand(tc.V,1);

            % a synthetic parcellation that does NOT cover the surface -- the real ones do
            % not either, and a denominator that assumes it does is silently wrong
            memb = false(tc.nROI, tc.V);
            for r = 1:tc.nROI, memb(r, (r-1)*12 + (1:12)) = true; end
            atlas = struct('Name','synthetic', 'nScout', tc.nROI, ...
                'Label', {compose('R%d', 1:tc.nROI)}, ...
                'Hemi',  {repmat({'L'},1,tc.nROI)}, ...
                'Membership', sparse(memb), 'nV', tc.V, ...
                'Vertices', {arrayfun(@(r) ((r-1)*12+(1:12))', 1:tc.nROI, 'UniformOutput', false)});

            tc.B = struct('pager', rheome.pagedrecording(f,'PageLength',400,'Overlap',100,'Precision','single'), ...
                'Kc', randn(tc.Ks, tc.C), 'Phi', Phi, 'Lambda', Lambda, 'fs', tc.fs, ...
                'gfb', rheome.graphfilterbank(max(Lambda),'NumFilters',3), 'wVert', wV, ...
                'atlas', atlas, 'name','synthetic', 'freqLimits',[10 120], ...
                'voicesOct', 6, 'chSel', (1:tc.C).', ...
                'surface', struct('Vertices',randn(tc.V,3),'Faces',[1 2 3],'nV',tc.V));
        end
    end

    methods (Test)

        function gramEnergyEqualsDirectPerFrameComputation(tc)
            % THE anchor. E_roi,g(win) must equal summing, frame by frame over the window,
            % the area-weighted |vorticity|^2 inside the ROI at that scale.
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            % ⚠ THE SAME METHOD as the extractor used. flowfeatures defaults to 'bandpass';
            % a reference built with 'cwt' has different coefficients, so the comparison
            % would fail for a reason that has nothing to do with the Gram algebra.
            fp  = rheome.flowpage(tc.B, 'Band', fe.BandLimits{1}, 'PageIndex', 1, ...
                           'FreqLimits', tc.B.freqLimits, 'VoicesPerOctave', tc.B.voicesOct, ...
                           'Method', fe.Method);

            w = tc.B.wVert;  Phi = tc.B.Phi;
            for iw = [1 2]
                idx = out.WindowFrames{iw, 1};
                for g = [1 3]
                    y = Phi * (double(fp.ScaleGains(:,g)) .* double(fp.Coefficients(:, idx)));
                    e = abs(y).^2;                                  % [V x nWinFrames]
                    for r = 1:tc.nROI
                        m = full(tc.B.atlas.Membership(r,:)).';
                        ref = sum(sum(e .* (w .* m), 1));
                        tc.verifyEqual(out.Energy(iw,1,g,r), ref, 'RelTol', 1e-6);
                    end
                end
            end
        end

        function windowsTileTheCoreAndNeverTheMargin(tc)
            % A window straddling the margin would mix this page's samples with the
            % neighbours' AND with cone-contaminated ones, and nothing downstream could tell.
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha'}, 'WindowSec', 0.25);
            out = extract(fe, 2);
            fp  = rheome.flowpage(tc.B, 'Band', fe.BandLimits{1}, 'PageIndex', 2, ...
                           'FreqLimits', tc.B.freqLimits, 'VoicesPerOctave', tc.B.voicesOct, ...
                           'Method', fe.Method);
            all_ = [];
            for iw = 1:out.NumWindows
                idx = out.WindowFrames{iw,1};
                tc.verifyTrue(all(fp.Core(idx)), 'a window reached into the margin');
                all_ = [all_, idx]; %#ok<AGROW>
            end
            tc.verifyEqual(numel(unique(all_)), numel(all_), 'windows overlap');
        end

        function everyWholeWindowInTheCoreIsUsed(tc)
            % ⚠ The core spans one sample MORE than (t_last - t_first): 4000 samples at
            % 400 Hz cover 9.9975 s between timestamps but 10.0 s of signal. Flooring the
            % timestamp span dropped a whole window per page -- 10% of the recording, with
            % nothing to indicate it.
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            fp  = rheome.flowpage(tc.B, 'Band', fe.BandLimits{1}, 'PageIndex', 1, ...
                'FreqLimits', tc.B.freqLimits, 'VoicesPerOctave', tc.B.voicesOct, ...
                'Method', fe.Method);
            coreDur = sum(fp.Core) / fp.Rate;
            tc.verifyEqual(out.NumWindows, floor(coreDur / fe.WindowSec + 1e-9));
            used = sum(out.WindowSamples(:,1));
            tc.verifyGreaterThan(used / sum(fp.Core), 0.98);   % ~all core samples used
        end

        function modeSpectrumIsTheDiagonalOfTheWindowCovariance(tc)
            % The full-resolution scale axis, free from the same Gwin.
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            tc.verifySize(out.ModeSpectrum, [out.NumWindows, 1, tc.Ks]);
            fp = rheome.flowpage(tc.B, 'Band', fe.BandLimits{1}, 'PageIndex', 1, ...
                'FreqLimits', tc.B.freqLimits, 'VoicesPerOctave', tc.B.voicesOct, ...
                'Method', fe.Method);
            idx = out.WindowFrames{2,1};
            ref = sum(abs(double(fp.Coefficients(:, idx))).^2, 2);
            tc.verifyEqual(squeeze(out.ModeSpectrum(2,1,:)), ref, 'RelTol', 1e-6);
            % and it must sum to the total
            tc.verifyEqual(sum(out.ModeSpectrum(2,1,:)), out.TotalEnergy(2,1), 'RelTol', 1e-9);
        end

        function shapeFollowsTheFourAxes(tc)
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha','beta'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            tc.verifySize(out.Energy, [out.NumWindows, 2, out.NumScales, tc.nROI]);
            tc.verifyEqual(out.Bands, {'alpha','beta'});
            tc.verifyEqual(numel(out.ROILabel), tc.nROI);
            tc.verifyEqual(numel(out.ScaleMM), out.NumScales);
        end

        function normalisersAreCarriedNotBakedIn(tc)
            % Raw energy is not comparable across subjects, but WHICH normalisation is right
            % depends on the question -- a table that has already divided cannot be undivided.
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            tc.verifySize(out.TotalEnergy,   [out.NumWindows, 1]);
            tc.verifySize(out.ScaleEnergy,   [out.NumWindows, 1, out.NumScales]);
            tc.verifySize(out.WindowSamples, [out.NumWindows, 1]);
            tc.verifyEqual(out.ROIArea(:), full(tc.B.atlas.Membership * tc.B.wVert), 'RelTol', 1e-12);
            % whole-cortex per-scale energy must be >= the sum over ROIs, since the
            % parcellation does not cover the surface
            tc.verifyGreaterThanOrEqual(squeeze(out.ScaleEnergy(1,1,:)) + 1e-9, ...
                                        squeeze(sum(out.Energy(1,1,:,:), 4)));
        end

        function windowSampleCountFollowsEachBandsOwnRate(tc)
            % ⚠ MULTIRATE: each band has its own derived rate, so a fixed-DURATION window
            % holds a band-dependent number of SAMPLES. Comparing summed energy across bands
            % without dividing by this compares sample counts, not physiology.
            %
            % (In this synthetic bundle fs = 400 and the 30-samples-per-cycle ceiling is
            % 13 Hz, so every band clamps to decim = 1 and the rates coincide. The invariant
            % worth pinning is not that they differ but that each window's count follows its
            % OWN band's rate -- which is what stays true at 2400 Hz where they do differ.)
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha','gamma'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            tc.verifyEqual(numel(out.Rate), 2);
            for b = 1:2
                tc.verifyEqual(out.WindowSamples(2,b), round(fe.WindowSec * out.Rate(b)), ...
                    'AbsTol', 1);
            end
            % and the count must equal the number of frames actually used
            tc.verifyEqual(out.WindowSamples(2,1), numel(out.WindowFrames{2,1}));
        end

        function longFormTableHasOneRowPerCell(tc)
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha','beta'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            T = totable(fe, out);
            tc.verifyEqual(height(T), out.NumWindows * 2 * out.NumScales * tc.nROI);
            tc.verifyTrue(all(ismember({'window','tStart','band','scaleMM','roi','hemi', ...
                'energy','density','totalEnergy','roiArea','nSamples'}, T.Properties.VariableNames)));
            tc.verifyClass(T.band, 'categorical');
            tc.verifyClass(T.roi,  'categorical');
        end

        function densityDividesOutSamplesAndArea(tc)
            % Energy summed over samples and area is not comparable between ROIs of different
            % size or bands of different rate; density is.
            fe = rheome.flowfeatures(tc.B, 'Bands', {'alpha'}, 'WindowSec', 0.25);
            out = extract(fe, 1);
            ref = out.Energy(2,1,2,3) / (out.WindowSamples(2,1) * out.ROIArea(3));
            tc.verifyEqual(out.Density(2,1,2,3), ref, 'RelTol', 1e-12);
        end

        function unknownBandNameErrors(tc)
            tc.verifyError(@() rheome.flowfeatures(tc.B, 'Bands', {'ultra'}), 'flowfeatures:band');
        end

    end
end

% Author: Diellor Basha, 2026
