function file = rawrecording(name, rawMatFile, opts)
% IMPORT.RAWRECORDING  Native-rate BST-BIN -> the pageable v7.3 store.
%
%   file = rheome.import.rawrecording(name, rawMatFile)
%   file = rheome.import.rawrecording(name, rawMatFile, opts)
%
% Writes +data/<name>/recording_native.mat with the sensor matrix as a TOP-LEVEL variable,
% exactly like rheome.import.recording, but read from the ORIGINAL acquisition rather than the
% resample cached in study.mat. On the reference CTF data that is 2400 Hz against 600 Hz.
%
% ⭐ WHY THE NATIVE RATE MATTERS. The multirate plan gives each filter 30 samples per cycle
% of its own frequency, which needs 30*f <= fs. At 600 Hz that ceiling is 20 Hz and 16 of
% the 60 filters in a [1 60] Hz bank fall short. At 2400 Hz the ceiling is 80 Hz and NONE
% do. Better still, multirate at 2400 Hz costs LESS per page (0.94 GB) than flat at 600 Hz
% (1.28 GB) -- full resolution is cheaper than what it replaces.
%
% ⚠ THE SSP PROJECTORS ARE APPLIED, and that is not optional. rheome.io.read.rawbst applies them
% as Brainstorm does; they carry the cardiac and blink artefact cleaning. Verified against
% Brainstorm's own import of the same file: WITHOUT them r = 0.80, WITH them r = 0.999993
% (min 0.999911 over 270 MEG channels).
%
% ⚠ STORED AS single -- and that is LOSSLESS here, unlike the study-cache path. The .bst
% holds float32, so single is BIT-EXACT and halves the store (1.73 GB per 600 s subject
% instead of 3.46 GB). rheome.import.recording stores double because its source is double.
%
% ⚠ WRITTEN IN CHUNKS through a matfile handle, so peak memory is one chunk, not the whole
% record. The default 60 s chunk is ~35 MB at 300 x 2400.
%
% INPUTS:
%   name        cached dataset name (the folder must already exist)
%   rawMatFile  path to the Brainstorm raw descriptor (data_0raw_*.mat)
%   opts .chunkSamples  samples per write (default 60 s worth)
%        .channelFile   override the auto-located channel file
%
% See also: rheome.io.read.rawbst, rheome.import.recording, rheome.pagedrecording
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    name  = char(name);
    dsdir = fullfile(rheome.load.root(), name);
    if ~exist(dsdir, 'dir')
        error('import:rawrecording:missing', ...
            'No dataset folder for ''%s''. Import the surface/study first.', name);
    end

    ropts = struct();
    if isfield(opts,'channelFile'), ropts.channelFile = opts.channelFile; end

    H = rheome.io.read.rawbst(rawMatFile, [], ropts);          % header only
    nCh = H.nCh;  nT = H.nT;  sfreq = H.sfreq;
    if ~isfield(opts,'chunkSamples') || isempty(opts.chunkSamples)
        opts.chunkSamples = max(1, round(60 * sfreq));
    end
    chunk = min(double(opts.chunkSamples), nT);

    fprintf('rheome.import.rawrecording[%s]: %d ch x %d samples @ %g Hz (%.0f s)\n', ...
        name, nCh, nT, sfreq, nT/sfreq);

    file = fullfile(dsdir, 'recording_native.mat');
    if exist(file, 'file'), delete(file); end           % matfile append would keep old data

    % Create the store, then fill F chunk by chunk through the matfile handle so peak
    % memory is one chunk. The header is written first so a partial file is still readable.
    Time        = (0:nT-1) / sfreq;                     %#ok<NASGU>
    ChannelFlag = H.ChannelFlag;                        %#ok<NASGU>
    ChannelName = H.ChannelName;                        %#ok<NASGU>
    ChannelType = H.ChannelType;                        %#ok<NASGU>
    Comment     = H.Comment;                            %#ok<NASGU>
    Events = []; nAvg = 1; nProjector = 0;              %#ok<NASGU>
    source = struct('raw', rawMatFile, 'bst', H.BstFile, 'native', true);  %#ok<NASGU>
    builtin('save', file, 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
        'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', 'nProjector', ...
        '-v7.3', '-nocompression');

    m = matfile(file, 'Writable', true);
    m.F = zeros(nCh, 0, 'single');                      % declare F, then grow by assignment
    nP = 0;  t0 = tic;
    for a = 1:chunk:nT
        b = min(a + chunk - 1, nT);
        R = rheome.io.read.rawbst(rawMatFile, [a b], ropts);
        m.F(1:nCh, a:b) = single(R.F);
        nP = R.nProjector;
        fprintf('  %6.1f%%  samples %d..%d\r', 100*b/nT, a, b);
    end
    m.nProjector = nP;
    d = dir(file);
    fprintf('  done in %.1f s | %d projector(s) applied | %s (%.0f MB)\n', ...
        toc(t0), nP, file, d.bytes/1e6);
end

% Author: Diellor Basha, 2026
