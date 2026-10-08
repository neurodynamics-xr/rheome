function info = importsubject(name, protoDir, opts)
% SCALE.IMPORTSUBJECT  One Brainstorm-protocol subject -> the data cache the analyses load from.
%
%   info = rheome.scale.importsubject('sub02', '/tmp/bst')
%   info = rheome.scale.importsubject(name, protoDir, DurationS=300, Rate=600, K=1000, Kc=400)
%
% Writes, under rheome.load.root()/<name>/ (set RHEOME_DATA to put it on node-local disk):
%   surface.mat   rheome.import.surface on tess_cortex_pial_low (the leadfield's own surface)
%   study.mat     chan, hm (unconstrained os_meg leadfield), ncov (noisecov_full), rec
%   noise.mat     the same-session noise run, same rate and projector treatment as rec
%   bases.mat     rheome.import.bases(name, K, Kc)  -- rheome.load.bases picks it up as the richest LBO cache
%   atlas.mat     rheome.import.atlas(name)
%
% ⭐ THE RECORDING IS READ FROM THE RAW .bst, NOT AN IMPORTED BLOCK. The protocol holds only the
% raw link (data_0raw_*), e.g. at 1200 Hz; rheome.io.read.rawbst reads it with the SSP projectors applied,
% exactly as Brainstorm would. It is then resampled to Rate (600 Hz by default) with resample(), and cropped to the first DurationS seconds.
%
% ⚠ CHOICES THAT COULD NOT BE AVOIDED:
%   * duration: every subject is cropped to the first DurationS seconds (300 s by default) so that
%     durations are equal across cohorts whose rest runs differ in length.
%   * empty room: each subject's own same-session noise run is used, at the same rate as the rest
%     recording.
%   * resampling: MATLAB resample (Kaiser FIR), not Brainstorm's; same rate, different anti-alias
%     filter from a study resampled inside Brainstorm.
%
% Returns info: the located files (rheome.scale.locate), durations, rates, and the kernel the protocol
% ships (recorded for provenance; the analyses build their own plain MNE, SnrFixed = 3, as the
% reports do).
%
% See also: rheome.scale.locate, rheome.scale.run, rheome.import.surface, rheome.import.bases, rheome.io.read.rawbst
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        protoDir (1,:) char
        opts.DurationS (1,1) double = 300
        opts.Rate (1,1) double = 600
        opts.K (1,1) double = 1000
        opts.Kc (1,1) double = 400
        opts.Sub char = ''
    end
    L = rheome.scale.locate(protoDir, opts.Sub);
    d = fullfile(rheome.load.root(), name);  if ~exist(d, 'dir'), mkdir(d); end
    info = struct('name', name, 'files', L, 'opts', opts);
    t0 = tic;

    rheome.import.surface(name, L.cortex);

    chan = rheome.io.read.channel(L.channel);                    %#ok<NASGU>
    hm   = rheome.io.read.headmodel(L.headmodel);
    ncov = rheome.io.read.noisecov(L.noisecov);                  %#ok<NASGU>
    [rec, info.rest] = i_raw(L.restRaw, opts);            %#ok<ASGLU>
    studyDir = L.restDir;  dataName = L.restRaw;           %#ok<NASGU>
    builtin('save', fullfile(d, 'study.mat'), 'chan', 'hm', 'rec', 'ncov', 'studyDir', 'dataName', '-v7.3');
    info.nChannels = rec.nCh;  info.leadfield = size(hm.Gain);
    clear rec hm

    info.noise = struct('file', '', 'durationS', 0);
    if ~isempty(L.noiseRaw)
        [nrec, info.noise] = i_raw(L.noiseRaw, opts);      %#ok<ASGLU>
        builtin('save', fullfile(d, 'noise.mat'), 'nrec', '-v7.3');
        clear nrec
    end

    rheome.import.bases(name, opts.K, opts.Kc);
    rheome.import.atlas(name);

    k = dir(fullfile(L.restDir, 'results_*KERNEL*.mat'));  k = k(~startsWith({k.name}, '._'));
    info.shippedKernels = string({k.name});
    info.importSeconds = toc(t0);
end

function [rec, meta] = i_raw(rawFile, opts)
% the raw link -> a study-style rec struct at opts.Rate, first opts.DurationS seconds
    H  = rheome.io.read.rawbst(rawFile);
    fs = H.sfreq;  nS = min(H.nT, round(opts.DurationS * fs));
    R  = rheome.io.read.rawbst(rawFile, [1 nS]);
    [p, q] = rat(opts.Rate / fs);
    F = resample(double(R.F)', p, q)';                     % column-wise over channels
    rec = struct('F', F, 'Time', (0:size(F,2)-1) / opts.Rate, 'ChannelFlag', R.ChannelFlag(:), ...
                 'nCh', size(F,1), 'nT', size(F,2), 'sfreq', opts.Rate, 'Comment', R.Comment, ...
                 'Events', [], 'nAvg', 1, 'ChannelName', {R.ChannelName}, 'ChannelType', {R.ChannelType});
    meta = struct('file', rawFile, 'nativeRate', fs, 'durationS', nS / fs, ...
                  'availableS', H.nT / fs, 'nProjector', R.nProjector);
end

% Author: Diellor Basha, 2026
