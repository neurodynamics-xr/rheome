function [E, cost] = bursts(db, opts)
% SELECT.BURSTS  Two-phase event query: prune on the stored envelope, refine on raw samples.
%
%   [E, cost] = rheome.select.bursts(db, Band=j, Threshold=theta, MinDuration=0.3, Channels=[], Margin=1.2, Halo=[])
%
% Phase 1 (atlas only): rheome.select.frames on envMax at the band's natural level with the
% threshold divided by Margin, keeping runs of MinDuration. The stored envMax is the max
% over the member's decimated instants (oversample 2), which can sit below the full-rate
% envelope's peak; Margin covers that gap so no event is lost (design §3.5; measured in
% tests/tSelectBursts.m against a brute-force scan).
% Phase 2 (raw): for each run, read the channel's samples over the run plus a halo of
% the band's support on each side, evaluate THE SAME MEMBERS of THE SAME BANK on that
% span at full rate (a timefilterbank at the span's length has the same members: the grid
% is anchored at an absolute frequency), take the envelope as the max over the band's
% members, and apply theta and MinDuration at sample resolution.
%
% OUTPUT E: recording_id, channel_id, band_id, t_on, t_off, duration, peak, t_peak, censored
%        (censored: within the band's support of the record edge)
%        cost: phase-1 cost, raw samples read, their fraction of a full scan of the queried
%        channels (spans are merged per channel, so it never exceeds 1), seconds
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        opts.Band        (1,1) double
        opts.Threshold   (1,1) double {mustBePositive}
        opts.MinDuration (1,1) double {mustBeNonnegative} = 0
        opts.Channels    double = []
        opts.Margin      (1,1) double {mustBeGreaterThanOrEqual(opts.Margin, 1)} = 1.2
        opts.Halo        double = []
    end
    t0 = tic;
    b = opts.Band;  bands = db.bands;  meta = db.meta;  fs = meta.fs;  nT = meta.nT;
    Lb = bands.naturalLevel(b);
    q = rheome.select.query(Stat="envMax", Op=">=", Threshold=opts.Threshold / opts.Margin, Level=Lb, Band=b, ...
                     Channels=opts.Channels, MinDuration=opts.MinDuration);
    [R, c1] = rheome.select.frames(db, q);
    halo = opts.Halo;  if isempty(halo), halo = 2 * bands.tSupport(b); end
    ext = db.grid.tExtent(Lb+1);
    E = table('Size', [0 9], 'VariableTypes', {'string','double','double','double','double','double','double','double','logical'}, ...
              'VariableNames', {'recording_id','channel_id','band_id','t_on','t_off','duration','peak','t_peak','censored'});
    samplesRead = 0;
    if isempty(R)
        cost = i_cost(c1, 0, nT, toc(t0));  return
    end
    if ~ismember('run_id', R.Properties.VariableNames)
        R.run_id = (1:height(R))';  R.run_start = R.k;  R.run_length = ones(height(R), 1);
    end
    runs = unique(R(:, {'channel_id','run_id','run_start','run_length'}), 'rows');
    V = meta.cfg.VoicesPerOctave;  a = meta.cfg.Anchor;
    chans = unique(runs.channel_id)';
    for ch = chans
        rc = runs(runs.channel_id == ch, :);
        tA = (rc.run_start - 1) * ext;  tB = (rc.run_start + rc.run_length - 1) * ext;
        s1 = max(1, floor((tA - halo) * fs) + 1);  s2 = min(nT, ceil((tB + halo) * fs));
        % merge overlapping spans, so no sample is read twice for a channel
        [s1, o] = sort(s1);  s2 = s2(o);  tA = tA(o);  tB = tB(o);
        spans = [s1(1) s2(1)];  cores = {[tA(1) tB(1)]};
        for i = 2:numel(s1)
            if s1(i) <= spans(end, 2) + 1
                spans(end, 2) = max(spans(end, 2), s2(i));  cores{end} = [cores{end}; tA(i) tB(i)];
            else
                spans(end+1, :) = [s1(i) s2(i)];  cores{end+1} = [tA(i) tB(i)]; %#ok<AGROW>
            end
        end
        for i = 1:size(spans, 1)
            x = i_read(db, ch, spans(i, 1), spans(i, 2));  samplesRead = samplesRead + numel(x);
            env = i_envelope(x, fs, V, a, bands, b);
            tt = (spans(i, 1) - 1 + (0:numel(x)-1)') / fs;
            core = false(numel(x), 1);
            for r = 1:size(cores{i}, 1)
                core = core | (tt >= cores{i}(r, 1) - halo/2 & tt < cores{i}(r, 2) + halo/2);   % not the span's own edges
            end
            on = env >= opts.Threshold & core;
            d = diff([0; on; 0]);  starts = find(d == 1);  stops = find(d == -1) - 1;
            for e = 1:numel(starts)
                dur = (stops(e) - starts(e) + 1) / fs;
                if dur < opts.MinDuration, continue; end
                [pk, ip] = max(env(starts(e):stops(e)));
                tOn = tt(starts(e));  tOff = tt(stops(e)) + 1/fs;
                cens = tOn < bands.tSupport(b) || tOff > nT/fs - bands.tSupport(b);
                E(end+1, :) = {string(db.recording_id), ch, b, tOn, tOff, dur, pk, tt(starts(e) + ip - 1), cens}; %#ok<AGROW>
            end
        end
    end
    E = unique(E, 'rows');                                      % a run's halo may see its neighbour's event
    E = sortrows(E, {'channel_id','t_on'});
    cost = i_cost(c1, samplesRead, nT * numel(chans), toc(t0));
end

function x = i_read(db, ch, s1, s2)
    pr = rheome.pagedrecording(db.meta.source, 'Channels', db.meta.iChannel(ch), 'Precision', 'double');
    x = read(pr, s1, s2).';
end

function env = i_envelope(x, fs, V, a, bands, b)
% the band's members of the anchored bank on this span, at full rate: max |w_m(t)|
    n = numel(x);
    tfb = rheome.timefilterbank(n, 'SamplingFrequency', fs, 'VoicesPerOctave', V, 'Anchor', a);
    [H, f] = freqz(tfb);  fc = centerFrequencies(tfb);  kd = kinds(tfb);
    switch bands.kind{b}
        case 'octave', mem = find(strcmp(kd, 'band') & fc >= bands.fLo(b) - 1e-9 & fc < bands.fHi(b) - 1e-9);
        case 'above',  mem = find(strcmp(kd, 'highpass'));
        case 'below',  mem = find(strcmp(kd, 'lowpass'));
    end
    X = fft(x);
    env = zeros(n, 1);
    for m = mem
        h = zeros(n, 1);  h(1:numel(f)) = H(m, :);
        env = max(env, abs(ifft(X .* h)));
    end
end

function cost = i_cost(c1, samplesRead, fullScan, secs)
% rawFraction: samples read over the samples a full scan of the queried channels would read
    cost = struct('phase1', c1, 'rawSamplesRead', samplesRead, 'rawFraction', samplesRead / fullScan, 'seconds', secs);
end
% Author: Diellor Basha, 2026
