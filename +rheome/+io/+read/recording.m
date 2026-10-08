function R = recording(DataFile)
% IO.READ.RECORDING  Read a Brainstorm data (recordings) file.
%
%   R = rheome.io.read.recording(DataFile)
%
% Loads a sensor recording (evoked average or continuous block).
%
% OUTPUT (struct R):
%   .F           [nCh x nT] sensor data
%   .Time        [1 x nT]  time vector (seconds)
%   .ChannelFlag [nCh x 1] +1 good / -1 bad
%   .nCh, .nT    counts
%   .sfreq       sampling frequency (Hz), derived from Time
%   .Comment, .Events, .nAvg
%
% Author: Diellor Basha, 2026

    raw = load(DataFile);
    if ~isfield(raw, 'F') || isempty(raw.F)
        error('io:read:recording:noF', 'No F matrix in %s', DataFile);
    end
    R = struct();
    R.Time        = i_get(raw, 'Time', []);
    R.ChannelFlag = i_get(raw, 'ChannelFlag', []);
    if isstruct(raw.F)
        % RAW link (continuous file): F is an sFile descriptor, not a data matrix. The binary
        % cannot be read from this pure-MATLAB module, so leave F empty and keep only the header
        % fields (ChannelFlag/Time). Data-independent uses -- e.g. building the inverse/flow
        % kernels -- are unaffected; streaming reads happen downstream (Brainstorm).
        warning('io:read:recording:raw', 'Raw link file (F is an sFile struct); F left empty: %s', DataFile);
        R.F   = [];
        R.nCh = numel(R.ChannelFlag);
        R.nT  = 0;
    else
        R.F   = double(raw.F);
        R.nCh = size(R.F, 1);
        R.nT  = size(R.F, 2);
    end
    if isempty(R.ChannelFlag), R.ChannelFlag = ones(R.nCh, 1); end
    R.sfreq       = [];
    if ~isempty(R.F) && numel(R.Time) > 1, R.sfreq = 1 / mean(diff(R.Time)); end
    R.Comment     = i_get(raw, 'Comment', '');
    R.Events      = i_get(raw, 'Events', []);
    R.nAvg        = i_get(raw, 'nAvg', 1);
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
