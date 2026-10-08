function E = spindle(X, fs, opts)
% DETECT.SPINDLE  Sleep spindles from the band envelope, with duration gating.
%
%   E = rheome.detect.spindle(X, fs)
%   E = rheome.detect.spindle(X, fs, opts)
%
% X is [nChan x nTime] (a row vector is one channel). Returns a TABLE, one row per event:
%   channel  startSample  endSample  startSec  durationSec  peakAmp  peakFreq
% Sample indices are in the ORIGINAL sampling grid, whatever rate the detection ran at.
%
% METHOD -- the conventional envelope detector. Band-pass to the spindle band, take the
% analytic envelope, smooth it, and mark excursions above a HIGH threshold, extended out to
% where they fall below a LOW one; then gate on duration and merge near neighbours. The two
% thresholds matter: a single one makes the measured duration a function of where the
% threshold sits, because the envelope's flanks are shallow.
%
% ⚠ THE THRESHOLD IS A ROBUST BASELINE, NOT A PERCENTILE, and that was a real bug. A
% percentile fails whenever events are SPARSE: in a minute of signal holding one 1 s spindle,
% the 90th percentile of the envelope is a NOISE level, so the detector fires continuously and
% reports plausible-looking events made entirely of background. Thresholds are therefore
% median + k * robust SD (MAD-based), which does not move when the event rate changes.
%
% ⚠ THE BASELINE IS GLOBAL, AND THIS DECIDES WHETHER STAGE CONTRAST SURVIVES. It is taken
% over the whole input (or opts.baseline), never per segment. Thresholding each segment
% against ITSELF finds "spindles" in wake at the same rate as in N2, because wake's own
% background is low -- the detector then measures its own threshold rather than the data.
% Pass the same recording, or an explicit baseline mask, for every stage you intend to compare.
%
% ⚠ THE BOUNDARY THRESHOLD GOVERNS THE FALSE-POSITIVE RATE, NOT THE DETECTION ONE -- which
% is not the intuitive way round. MEASURED on 5 minutes of narrowband noise plus 21 planted
% spindles at 2.5x background:
%
%     hiSD  loSD | hits/21 | false per min on noise
%      3.0   1.5 |   21    |   0.40
%      3.0   2.5 |   21    |   0.00
%      6.0   1.5 |   21    |   0.00 ... but 5.0/1.5 still gave 0.20
%
% Raising hiSD barely helps because a noise excursion that touches it is still EXTENDED down
% to loSD, and it is that extension which carries the run past the duration gate. Lowering
% loSD lengthens noise events into acceptance. Sensitivity was 21/21 at every setting tried,
% so loSD is the parameter with a cost, and 2.5 is the default for that reason.
%
% ⚠ REPORTED DURATION IS AMPLITUDE-DEPENDENT near the gate. Smoothing the envelope and
% extending each event down to the boundary threshold WIDENS it, and the wider the louder:
% measured, a 0.2 s burst at 109 robust-SD above baseline stays above the boundary for 2.4 s
% and is reported as a 0.80 s event. At realistic spindle amplitudes (~2.5x background) the
% gate behaves; do not read a duration off a very high-amplitude event.
%
% ⚠ IT WILL FIRE ON ARTEFACT AND LOOK FINE DOING IT. Eye movement and muscle produce
% envelope excursions of spindle-like duration. The check that catches it is TOPOGRAPHY, not
% the event list: real spindles are central-parietal, blinks are frontopolar. Always look at
% where the detections are before believing how many there are.
%
% ⚠ DECIMATES INTERNALLY. A 11-16 Hz band at 1000 Hz sits at 2-3% of Nyquist where a
% low-order IIR is poorly conditioned, and filtering 28 million samples per channel is slow
% for no benefit. The signal is resampled to opts.workRate (default 100 Hz, ample for a
% 16 Hz band) and every reported index is mapped back to the original grid.
%
% INPUTS
%   opts .band       [11 16] Hz          .minDur   0.5 s     .maxDur  3.0 s
%        .hiSD       3.0 detection threshold, median + k*robustSD
%        .loSD       2.5 boundary threshold, same units
%        .mergeGap   0.1 s  merge events closer than this
%        .smoothSec  0.2 s  envelope smoothing
%        .workRate   100 Hz internal detection rate
%        .baseline   logical [1 x nTime] samples used for the percentiles (default: all)
%        .envelope   false; true = X IS ALREADY AN ENVELOPE (e.g. a source-space |J^| pooled over
%                    a tile): the band-pass and Hilbert are skipped, everything after them --
%                    smoothing, the two robust thresholds, extension, merging, gating -- is the
%                    same, and peakFreq is NaN (there is no carrier left to read it from).
%                    ⚠ .band is then only a label; the envelope's own band was set upstream.
%
% See also: rheome.detect.vortex, rheome.detect.track, rheome.filters.firbandpass
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    d = @(f,v) i_def(opts,f,v);
    band     = d('band',[11 16]);   minDur = d('minDur',0.5);  maxDur = d('maxDur',3.0);
    hiSD     = d('hiSD',3.0);       loSD   = d('loSD',2.5);
    mergeGap = d('mergeGap',0.1);   smoothSec = d('smoothSec',0.2);
    workRate = d('workRate',100);   baseline  = d('baseline',[]);
    isEnv    = d('envelope',false);

    if isvector(X), X = X(:).'; end
    [nCh, nT] = size(X);
    if ~isEnv && band(2) >= 0.45*fs
        error('detect:spindle:band', ...
            'Band [%g %g] Hz is not representable at %g Hz (need band(2) < %.1f Hz).', ...
            band(1), band(2), fs, 0.45*fs);
    end

    % ---- decimate to the working rate ----
    wr = min(workRate, fs);
    [p,q] = rat(wr/fs);
    Xw = resample(double(X).', p, q).';
    fsw = fs * p/q;
    nW  = size(Xw,2);
    bw = [];
    if ~isempty(baseline)
        bw = resample(double(baseline(:)), p, q) > 0.5;
        bw = bw(1:min(numel(bw),nW));  bw(end+1:nW) = false;
    end

    % ---- band envelope ----
    if isEnv
        env = max(Xw, 0);  Yb = [];                         % resampling can ring below zero
    else
        dfl = designfilt('bandpassiir','FilterOrder',6, ...
            'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',min(band(2),0.45*fsw), ...
            'SampleRate',fsw);
        if ~isstable(dfl)
            error('detect:spindle:filter','No stable filter for [%g %g] Hz at %g Hz.', ...
                band(1), band(2), fsw);
        end
        Yb  = filtfilt(dfl, Xw.').';
        env = abs(hilbert(Yb.')).';
    end
    k   = max(1, round(smoothSec*fsw));
    env = movmean(env, k, 2);

    rows = {};
    for c = 1:nCh
        e = env(c,:);
        ref = e;  if ~isempty(bw), ref = e(bw); end
        if isempty(ref) || all(ref==0), continue; end
        % robust centre and spread: the median and MAD are unmoved by the events themselves,
        % where a mean or a percentile would be dragged up by them
        med = median(ref);
        rsd = 1.4826 * median(abs(ref - med));
        if ~(rsd > 0), continue; end
        hi = med + hiSD*rsd;  lo = med + loSD*rsd;

        above = e >= lo;
        [s0, s1] = i_runs(above);
        keep = arrayfun(@(a,b) any(e(a:b) >= hi), s0, s1);   % must reach the HIGH threshold
        s0 = s0(keep); s1 = s1(keep);
        if isempty(s0), continue; end

        % merge near neighbours before gating -- a spindle briefly dipping under the low
        % threshold would otherwise be split into two short events and then both rejected
        g = round(mergeGap*fsw);
        m0 = s0(1); m1 = s1(1); M0=[]; M1=[];
        for i = 2:numel(s0)
            if s0(i) - m1 <= g, m1 = s1(i);
            else, M0(end+1)=m0; M1(end+1)=m1; m0=s0(i); m1=s1(i); end %#ok<AGROW>
        end
        M0(end+1)=m0; M1(end+1)=m1;

        dur = (M1-M0+1)/fsw;
        ok  = dur >= minDur & dur <= maxDur;
        M0=M0(ok); M1=M1(ok); dur=dur(ok);
        for i = 1:numel(M0)
            [pk, ipk] = max(e(M0(i):M1(i)));
            f = NaN;  if ~isEnv, f = i_peakfreq(Yb(c, M0(i):M1(i)), fsw, band); end
            % map back to the ORIGINAL grid
            a = round((M0(i)-1)*fs/fsw) + 1;
            b = min(round((M1(i)-1)*fs/fsw) + 1, nT);
            rows(end+1,:) = {c, a, b, (a-1)/fs, dur(i), pk, f, ...
                             (M0(i)+ipk-2)/fsw}; %#ok<AGROW>
        end
    end

    if isempty(rows)
        E = table('Size',[0 8], ...
            'VariableTypes',{'double','double','double','double','double','double','double','double'}, ...
            'VariableNames',{'channel','startSample','endSample','startSec','durationSec', ...
                             'peakAmp','peakFreq','peakSec'});
        return;
    end
    E = cell2table(rows, 'VariableNames', {'channel','startSample','endSample','startSec', ...
        'durationSec','peakAmp','peakFreq','peakSec'});
    E = sortrows(E, {'startSec','channel'});
end

% ---- runs of true ----
function [s0, s1] = i_runs(m)
    m = m(:).';  d = diff([false m false]);
    s0 = find(d==1);  s1 = find(d==-1) - 1;
end

% ---- dominant frequency inside the band, from the segment's own spectrum ----
function f = i_peakfreq(seg, fs, band)
    n = numel(seg);
    if n < 8, f = NaN; return; end
    nf = 2^nextpow2(max(n, round(4*fs)));
    P = abs(fft(seg(:).' .* hann(n).', nf)).^2;
    fr = (0:nf-1)*fs/nf;
    m = fr >= band(1) & fr <= band(2);
    [~,i] = max(P(m));  fm = fr(m);  f = fm(i);
end

function v = i_def(o, f, v)
    if isfield(o,f) && ~isempty(o.(f)), v = o.(f); end
end

% Author: Diellor Basha, 2026
