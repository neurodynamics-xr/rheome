function R = rawbst(rawMatFile, samples, opts)
% IO.READ.RAWBST  Read a range of samples from a Brainstorm BST-BIN raw file.
%
%   R = rheome.io.read.rawbst(rawMatFile)                 % header only, no data read
%   R = rheome.io.read.rawbst(rawMatFile, [a b])          % samples a..b, projectors applied
%   R = rheome.io.read.rawbst(rawMatFile, [a b], opts)
%
% rawMatFile is the sFile DESCRIPTOR (data_0raw_*.mat) that sits beside the .bst binary.
% This gives the native acquisition rate -- 2400 Hz on the reference CTF data, against the 600 Hz
% resample cached in study.mat -- without any Brainstorm dependency.
%
% ⚠ THE FORMAT, ESTABLISHED BY MEASUREMENT (verified r = 1.0000, min 0.9999, against
% Brainstorm's own import of the same file; every one of these was wrong at least once):
%
%   * float32, LITTLE-ENDIAN (sFile.byteorder = 'l'). Big-endian returns NaN.
%   * values are ALREADY in Tesla and ALREADY CTF-compensated (currCtfComp = destCtfComp
%     = 3), so no gain calibration is needed. MEG std ~ 2.2e-13.
%   * ⭐ EPOCH-BLOCKED in header.epochsize blocks, CHANNEL-MAJOR inside each block --
%     and this holds EVEN THOUGH sFile.epochs IS EMPTY. A continuous file is still written
%     in blocks. Reading it as one flat channel-major run (which the offset formula in
%     in_fread_bst reads like, for a single epoch) gives correct-looking per-channel
%     MAGNITUDES and a correlation of ZERO. That failure mode is silent.
%
% ⚠ THE SSP/ICA PROJECTORS ARE PART OF THE DATA, NOT AN OPTION. Brainstorm applies
% `F = Projector * F` on every read (in_fread). These carry the artefact cleaning -- on
% a reference subject, cardiac and blink SSP components. Skipping them is not "raw data": it is
% measurably dirtier data. MEASURED against Brainstorm's import of the same samples:
%
%     epoch-blocked, no projector ....... r = 0.8028
%     epoch-blocked, WITH projector ..... r = 1.0000   (min 0.9999, 100% of ch > 0.99)
%
% Only projectors with Status == 1 are applied (0 = disabled, 2 = already applied to the
% file). Applying a disabled one attenuates real signal.
%
% INPUTS:
%   rawMatFile  path to data_0raw_*.mat (the sFile descriptor)
%   samples     [a b] 1-based inclusive sample range; omit for a header-only read
%   opts .projector  true (default) | false  -- false is for VALIDATION ONLY
%        .channelFile  path to the channel file (default: channel_*.mat beside the raw)
%        .precision    'double' (default) | 'single'
%
% OUTPUT (struct R):
%   .F [nCh x n] (empty for a header-only read)   .Time [1 x n]   .samples [a b]
%   .nCh .nT .sfreq .ChannelFlag .ChannelName .ChannelType .Comment
%   .nProjector  active projectors applied        .BstFile  the binary actually read
%
% See also: rheome.import.rawrecording, rheome.pagedrecording, rheome.io.read.recording
%
% Author: Diellor Basha, 2026

    if nargin < 2, samples = []; end
    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'projector') || isempty(opts.projector), opts.projector = true;  end
    if ~isfield(opts,'precision') || isempty(opts.precision), opts.precision = 'double'; end
    if ~exist(rawMatFile,'file')
        error('io:read:rawbst:notFound', 'No such raw descriptor: %s', rawMatFile);
    end

    D = builtin('load', rawMatFile);
    if ~isfield(D,'F') || ~isstruct(D.F)
        error('io:read:rawbst:notRaw', ...
            '%s is not a raw LINK file (F is not an sFile struct).', rawMatFile);
    end
    sF = D.F;  h = sF.header;
    nCh = double(h.nchannels);  nS = double(h.nsamples);  ep = double(h.epochsize);
    hdr = double(h.hdrsize);

    R = struct();
    R.nCh = nCh;  R.nT = nS;
    R.sfreq = double(sF.prop.sfreq);
    R.ChannelFlag = i_get(sF, 'channelflag', ones(nCh,1));
    R.ChannelFlag = R.ChannelFlag(:);
    R.Comment  = i_get(D, 'Comment', '');
    R.samples  = [];  R.Time = [];  R.F = [];  R.nProjector = 0;

    % the binary: trust the descriptor's path, fall back to a sibling .bst
    bst = i_get(sF, 'filename', '');
    if isempty(bst) || ~exist(bst,'file')
        g = dir(fullfile(fileparts(rawMatFile), '*.bst'));
        g = g(~startsWith({g.name}, '._'));
        if isempty(g)
            error('io:read:rawbst:noBinary', 'No .bst binary found beside %s', rawMatFile);
        end
        bst = fullfile(g(1).folder, g(1).name);
    end
    R.BstFile = bst;

    % channel metadata + projectors
    if ~isfield(opts,'channelFile') || isempty(opts.channelFile)
        g = dir(fullfile(fileparts(rawMatFile), 'channel_*.mat'));
        g = g(~startsWith({g.name}, '._'));
        if isempty(g), opts.channelFile = ''; else, opts.channelFile = fullfile(g(1).folder, g(1).name); end
    end
    C = struct('Channel', [], 'Projector', [], 'MegRefCoef', []);
    if ~isempty(opts.channelFile) && exist(opts.channelFile,'file')
        C = builtin('load', opts.channelFile);
    end
    R.ChannelName = {};  R.ChannelType = {};
    if isfield(C,'Channel') && ~isempty(C.Channel)
        R.ChannelName = {C.Channel.Name};  R.ChannelType = {C.Channel.Type};
    end

    if isempty(samples), return; end                 % header-only read

    a = double(samples(1));  b = double(samples(end));
    if a < 1 || b > nS || b < a
        error('io:read:rawbst:range', ...
            'Sample range [%d %d] is outside 1..%d.', a, b, nS);
    end

    % ---- read whole epoch blocks, then slice ----
    % Each block is [nCh x ep] stored CHANNEL-MAJOR, so a block reshapes as [ep x nCh].'.
    e0 = floor((a-1)/ep);  e1 = floor((b-1)/ep);
    order = 'ieee-le';
    if isfield(sF,'byteorder') && strcmpi(sF.byteorder,'b'), order = 'ieee-be'; end
    fid = fopen(bst, 'r', order);
    if fid < 0, error('io:read:rawbst:open', 'Cannot open %s', bst); end
    cleanup = onCleanup(@() fclose(fid));

    blocks = cell(1, e1-e0+1);
    for e = e0:e1
        nThis = min(ep, nS - e*ep);                  % the final block may be short
        if fseek(fid, hdr + e*nCh*ep*4, 'bof') ~= 0
            error('io:read:rawbst:seek', 'Seek failed at epoch %d in %s', e, bst);
        end
        raw = fread(fid, [1 nCh*nThis], 'float32=>double');
        if numel(raw) < nCh*nThis
            error('io:read:rawbst:short', ...
                'Short read at epoch %d (%d of %d values).', e, numel(raw), nCh*nThis);
        end
        blocks{e-e0+1} = reshape(raw, nThis, nCh).';
    end
    B = [blocks{:}];
    s0 = a - e0*ep;
    F = B(:, s0 : s0 + (b-a));

    % ---- CTF gradient compensation, only if the file is not already at the target ----
    % Mirrors in_fread: F(iMeg,:) = F(iMeg,:) - MegRefCoef * F(iRef,:). On these files
    % currCtfComp == destCtfComp == 3, so this is a no-op -- kept so a differently
    % compensated file does not read silently wrong.
    if isfield(C,'MegRefCoef') && ~isempty(C.MegRefCoef) && ...
            ~isequal(sF.prop.currCtfComp, sF.prop.destCtfComp) && ~isempty(R.ChannelType)
        iMeg = find(strcmp(R.ChannelType, 'MEG'));
        iRef = find(strcmp(R.ChannelType, 'MEG REF'));
        F(iMeg,:) = F(iMeg,:) - C.MegRefCoef * F(iRef,:);
    end

    % ---- SSP / ICA projectors ----
    if opts.projector
        [P, nP] = i_buildprojector(i_get(C,'Projector',[]), nCh);
        if nP > 0, F = P * F; end
        R.nProjector = nP;
    end

    if strcmpi(opts.precision,'single'), F = single(F); end
    R.F       = F;
    R.samples = [a b];
    R.Time    = (a-1 : b-1) / R.sfreq;
end

% ---- Brainstorm's process_ssp2('BuildProjector', P, 1): product of (I - U*U') ----
function [P, nP] = i_buildprojector(Projector, nCh)
    P = [];  nP = 0;
    if isempty(Projector), return; end
    P = eye(nCh);
    for k = 1:numel(Projector)
        pk = Projector(k);
        if ~isfield(pk,'Components') || isempty(pk.Components), continue; end
        if ~isfield(pk,'Status') || pk.Status ~= 1, continue; end     % 0 off, 2 already applied
        mask = true(1, size(pk.Components,2));
        if isfield(pk,'CompMask') && ~isempty(pk.CompMask), mask = pk.CompMask ~= 0; end
        U = pk.Components(:, mask);
        if isempty(U), continue; end
        P  = (eye(nCh) - U*U') * P;
        nP = nP + 1;
    end
end

function v = i_get(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
