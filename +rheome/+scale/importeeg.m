function info = importeeg(name, protoDir, edfFile, eventsFile, opts)
% SCALE.IMPORTEEG  A sleep EEG subject -> the data cache the analyses load from: template forward model + scored NREM.
%
%   info = rheome.scale.importeeg('anphysleep_EPCTL03', '/tmp/bst', edf, eventsTsv)
%   info = rheome.scale.importeeg(name, protoDir, edf, eventsTsv, Rate=100, Stages=["N2" "N3"], BadFrac=0.5)
%
% The EEG counterpart of rheome.scale.importsubject, for the MS1 positive control (G15; nsp cf-slowosc).
% protoDir is a Brainstorm protocol holding the EEG FORWARD MODEL (Data D7, b3513816): a channel file whose
% EEG channels carry the electrode positions, and an UNCONSTRAINED surface head model (OpenMEEG BEM). The
% recording itself is NOT read from the protocol but from the EEG-BIDS EDF and its events.tsv
% (rheome.io.read.sleepeeg), so the protocol needs no raw link.
%
% ⭐ ANPHYSLEEP HAS NO INDIVIDUAL MRI AND NO INDIVIDUAL DIGITISATION (Data 4e19b265 inventory): the anatomy is
%   the ICBM152 template and the electrodes the release's co-registered average (.pos, ALS, metres), so the
%   forward model is the SAME for every participant. info.forward records which head model was used
%   (its Comment and file), so the template choice is stated per subject in provenance.json.
% ⚠ The cortex is the head model's own SurfaceFile, resolved under protoDir/anat: a leadfield on another surface
%   would mis-index every vertex.
% ⚠ Bad channels: a channel the release's artefact matrix marks in more than BadFrac of the kept epochs is
%   flagged -1 for the whole night; the per-epoch marks stay in rec.epochs.artifact for the analysis.
% ⚠ Noise covariance: identity (sleep EEG has no noise recording; the reference choice is the analysis's: the
%   EEG rows are average-referenced in rheome.scale.sensors with Modality "EEG").
%
% Writes, under rheome.load.root()/<name>/: surface.mat, study.mat (chan, hm, rec, ncov), bases.mat, atlas.mat
% (when the cortex carries one). Returns info (files, forward, channels, bad channels, epochs, seconds).
%
% See also: rheome.io.read.sleepeeg, rheome.scale.importsubject, rheome.scale.measure_slowosc
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        protoDir (1,:) char
        edfFile (1,:) char
        eventsFile (1,:) char
        opts.Rate (1,1) double = 100
        opts.Stages string = ["N2" "N3"]
        opts.BadFrac (1,1) double = 0.5
        opts.K (1,1) double = 1000
        opts.Kc (1,1) double = 400
    end
    t0 = tic;
    hmf = dir(fullfile(protoDir, 'data', '**', 'headmodel_surf_*.mat'));  hmf = hmf(~startsWith({hmf.name}, '._'));
    if numel(hmf) ~= 1, error('scale:importeeg:headmodel', 'expected one surface head model under %s/data, found %d', protoDir, numel(hmf)); end
    hmFile = fullfile(hmf.folder, hmf.name);
    chf = dir(fullfile(hmf.folder, 'channel*.mat'));  chf = chf(~startsWith({chf.name}, '._'));
    if numel(chf) ~= 1, error('scale:importeeg:channel', 'expected one channel file next to %s', hmFile); end
    raw = load(hmFile, 'SurfaceFile', 'Comment');
    cortex = fullfile(protoDir, 'anat', char(raw.SurfaceFile));
    if ~isfile(cortex), error('scale:importeeg:surface', 'the head model''s surface %s is not in the protocol', raw.SurfaceFile); end

    d = fullfile(rheome.load.root(), name);  if ~exist(d, 'dir'), mkdir(d); end
    rheome.import.surface(name, cortex);
    chan = rheome.io.read.channel(fullfile(chf.folder, chf.name));
    hm = rheome.io.read.headmodel(hmFile);
    isE = strcmpi(chan.Type, 'EEG');
    hdr = edfinfo(edfFile);  inEdf = ismember(chan.Name, cellstr(string(hdr.SignalLabels)));
    keep = find(isE & inEdf);
    if isempty(keep), error('scale:importeeg:channels', 'no EEG channel of %s is in %s', chf.name, edfFile); end
    chan.Channel = chan.Channel(keep);  chan.Type = chan.Type(keep);  chan.Name = chan.Name(keep);  chan.nCh = numel(keep);
    hm.Gain = hm.Gain(keep, :);  hm.nCh = numel(keep);

    rec = rheome.io.read.sleepeeg(edfFile, eventsFile, Stages=opts.Stages, Rate=opts.Rate, Channels=string(chan.Name));
    bad = i_badfrac(rec.epochs.artifact, string(chan.Name));
    rec.ChannelFlag(bad > opts.BadFrac) = -1;
    ncov = struct('NoiseCov', eye(chan.nCh), 'FourthMoment', [], 'nSamples', []);
    studyDir = hmf.folder;  dataName = edfFile; %#ok<NASGU>
    builtin('save', fullfile(d, 'study.mat'), 'chan', 'hm', 'rec', 'ncov', 'studyDir', 'dataName', '-v7.3');

    rheome.import.bases(name, opts.K, opts.Kc);
    try, rheome.import.atlas(name); catch e, info.atlasError = string(e.message); end %#ok<NOSEMI>
    info.name = name;  info.files = struct('headmodel', hmFile, 'channel', fullfile(chf.folder, chf.name), 'cortex', cortex, ...
        'edf', edfFile, 'events', eventsFile);
    info.forward = struct('comment', string(raw.Comment), 'surface', string(raw.SurfaceFile), 'template', contains(lower(string(raw.SurfaceFile)), ["default" "icbm"]));
    info.nChannels = chan.nCh;  info.badChannels = string(chan.Name(bad > opts.BadFrac));
    info.nEpochs = height(rec.epochs);  info.stages = opts.Stages;  info.keptS = size(rec.F, 2) / rec.sfreq;
    info.importSeconds = toc(t0);
end

function f = i_badfrac(art, names)
% fraction of the kept epochs in which each channel is marked (n/a counts as unmarked)
    f = zeros(numel(names), 1);  n = numel(art);
    for k = 1:n
        if art(k) == "none" || art(k) == "n/a", continue, end
        f = f + ismember(names(:), strtrim(split(art(k), ",")));
    end
    f = f / max(n, 1);
end

% Author: Diellor Basha, 2026
