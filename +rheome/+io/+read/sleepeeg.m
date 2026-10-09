function rec = sleepeeg(edfFile, eventsFile, opts)
% IO.READ.SLEEPEEG  One overnight EEG-BIDS recording -> the scored NREM epochs, resampled, with their stages and artefact marks.
%
%   rec = rheome.io.read.sleepeeg(edf, eventsTsv)
%   rec = rheome.io.read.sleepeeg(edf, eventsTsv, Stages=["N2" "N3"], Rate=100, Channels=names, ChunkS=300)
%
% The reader of the MS1 positive control (G15, nsp cf-slowosc). AnphySleep's EEG-BIDS (Data 4e19b265) holds
% one ~7 h EDF+ at 1000 Hz per subject and an events.tsv with one row per scored 30-s epoch,
% trial_type "stage/<label>", and the release's artefact matrix as artifact_channels ("none", or the
% comma-separated EEG channels marked bad in that epoch). Only the epochs whose stage is in Stages are read:
% edfread's SelectedDataRecords takes them a chunk at a time, so the whole night is never in memory.
%
% ⭐ ONLY THE SCORED EPOCHS IN Stages ARE KEPT, CONCATENATED. rec.seg says where each contiguous piece of
%   the night sits in rec.F; a window that crosses a segment boundary crosses a gap in time and must be
%   dropped (rheome.scale.measure_slowosc does).
% ⚠ Each run of consecutive epochs is read and resampled in chunks of ChunkS seconds, and each chunk is its
%   own segment: resample()'s edge transient (a few samples) never lands inside a window that is kept.
% ⚠ Values are the EDF's physical values (µV for AnphySleep, channels.tsv); nothing is re-referenced here.
%
% OUTPUT rec: .F single [nCh x nS]  .sfreq  .ChannelName {nCh x 1}  .ChannelFlag [nCh x 1] (all 1)
%   .epochs table: onsetS (in the night), stage, artifact (string), first (column of rec.F), n (samples)
%   .seg [nSeg x 2] first and last column of each contiguous segment   .sourceRate  .file
%
% See also: rheome.scale.importeeg, rheome.scale.measure_slowosc
%
% Author: Diellor Basha, 2026

    arguments
        edfFile (1,:) char
        eventsFile (1,:) char
        opts.Stages string = ["N2" "N3"]
        opts.Rate (1,1) double {mustBePositive} = 100
        opts.Channels string = string.empty
        opts.ChunkS (1,1) double {mustBePositive} = 300
    end
    ev = readtable(eventsFile, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
    isEp = startsWith(ev.trial_type, "stage/");
    ep = ev(isEp, :);  ep.stage = extractAfter(ep.trial_type, "stage/");
    ep = ep(ismember(ep.stage, opts.Stages), :);
    if ~ismember('artifact_channels', ep.Properties.VariableNames), ep.artifact_channels = repmat("none", height(ep), 1); end
    ep.artifact_channels(ismissing(ep.artifact_channels)) = "n/a";

    info = edfinfo(edfFile);  recS = seconds(info.DataRecordDuration);
    lab = string(info.SignalLabels);
    if isempty(opts.Channels), sel = lab; else, sel = opts.Channels(:)'; end
    missing = setdiff(sel, lab);
    if ~isempty(missing), error('io:read:sleepeeg:channels', 'not in %s: %s', edfFile, strjoin(missing, ', ')); end
    fs0 = info.NumSamples(find(lab == sel(1), 1)) / recS;
    [p, q] = rat(opts.Rate / fs0);

    % runs of consecutive epochs, then chunks of at most ChunkS seconds
    on = ep.onset;  du = ep.duration;
    brk = [true; abs(on(2:end) - (on(1:end-1) + du(1:end-1))) > 1e-6];
    run = cumsum(brk);  F = cell(0, 1);  seg = zeros(0, 2);  first = zeros(height(ep), 1);  n = first;  col = 0;
    for r = 1:max([run; 0])
        ie = find(run == r);  t0 = on(ie(1));  t1 = on(ie(end)) + du(ie(end));
        for c0 = t0:opts.ChunkS:t1 - 1e-9
            c1 = min(c0 + opts.ChunkS, t1);
            recs = (floor(c0 / recS) + 1):ceil(c1 / recS);
            D = edfread(edfFile, 'SelectedSignals', sel, 'SelectedDataRecords', recs);
            X = zeros(numel(sel), numel(vertcat(D{:, 1}{:})));
            for k = 1:numel(sel), X(k, :) = vertcat(D{:, k}{:}); end
            i0 = round((c0 - (recs(1) - 1) * recS) * fs0) + 1;  i1 = round((c1 - (recs(1) - 1) * recS) * fs0);
            Y = resample(X(:, i0:i1)', p, q)';
            F{end+1, 1} = single(Y); %#ok<AGROW>
            seg(end+1, :) = col + [1 size(Y, 2)]; %#ok<AGROW>
            inC = ie(on(ie) >= c0 - 1e-9 & on(ie) < c1 - 1e-9);
            first(inC) = col + round((on(inC) - c0) * opts.Rate) + 1;  n(inC) = round(du(inC) * opts.Rate);
            col = col + size(Y, 2);
        end
    end
    rec = struct('F', [F{:}], 'sfreq', opts.Rate, 'ChannelName', {cellstr(sel(:))}, 'ChannelFlag', ones(numel(sel), 1), ...
        'epochs', table(on, ep.stage, ep.artifact_channels, first, n, 'VariableNames', {'onsetS','stage','artifact','first','n'}), ...
        'seg', seg, 'sourceRate', fs0, 'file', edfFile);
end

% Author: Diellor Basha, 2026
