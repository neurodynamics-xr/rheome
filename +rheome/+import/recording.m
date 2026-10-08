function file = recording(name)
% IMPORT.RECORDING  Surface a cached recording as a v7.3, chunk-readable store.
%
%   file = rheome.import.recording(name)
%
% Reads +data/<name>/study.mat and writes +data/<name>/recording.mat with the sensor
% matrix as a TOP-LEVEL variable, so a page can be read without loading the record:
%
%   m = matfile(file);   X = m.F(:, a:b);        % a genuine partial read
%
% ⚠ WHY THIS EXISTS AT ALL. study.mat already holds the same numbers, and is already
% v7.3 -- but as rec.F, NESTED inside a struct. MatFile supports '()' indexing only on
% TOP-LEVEL variables: `m.rec.F(:, a:b)` errors outright, and `m.rec` pulls the entire
% 977 MB struct into memory. Nesting therefore defeats paging completely. Surfacing F as
% its own variable is the whole content of this function.
%
% ⚠ WRITTEN -NOCOMPRESSION, DELIBERATELY. Measured on a 300 x 360000 double:
%
%     write compressed    14.6 s   833 MB     4 s page read  13 ms   ( 443 MB/s)
%     write nocompression  0.9 s   865 MB     4 s page read   2 ms   (2649 MB/s)
%
% 6-12x faster to read and 16x faster to write for 4% more disk. A page store is read far
% more often than it is written, so the trade is not close.
%
% ⚠ STORED IN DOUBLE, read in single. The store is a faithful cache of the recording, not
% of any one analysis; pagedrecording casts on read (default 'single'), which is where the
% precision choice belongs. Costs ~865 MB per 600 s subject, ON TOP of study.mat -- the
% cache is regenerable, so delete study.mat's copy only if disk actually bites.
%
% See also: rheome.pagedrecording, rheome.import.study, rheome.load.root
%
% Author: Diellor Basha, 2026

    name  = char(name);
    dsdir = fullfile(rheome.load.root(), name);
    sfile = fullfile(dsdir, 'study.mat');
    if ~exist(sfile, 'file')
        error('import:recording:missing', ...
            'No study cached for ''%s''. Run  rheome.import.study(''%s'', studyDir, dataName)  first.', ...
            name, name);
    end

    fprintf('rheome.import.recording[%s]: reading cached study\n', name);
    S   = builtin('load', sfile, 'rec', 'chan', 'studyDir', 'dataName');
    rec = S.rec;

    if ~isfield(rec, 'F') || isempty(rec.F)
        % rheome.io.read.recording leaves F empty for a raw LINK file (F is an sFile descriptor,
        % not a matrix). Catch it here rather than at the first page(), far from the cause.
        error('import:recording:noData', ...
            ['Cached recording for ''%s'' carries no data matrix (raw link file). ' ...
             'Import the continuous block into its own Brainstorm condition and re-run ' ...
             'rheome.import.study with that data file.'], name);
    end

    F     = double(rec.F);
    nCh   = size(F, 1);
    nT    = size(F, 2);
    Time  = i_get(rec, 'Time', (0:nT-1));
    Time  = double(Time(:)).';
    sfreq = i_get(rec, 'sfreq', []);
    if isempty(sfreq) && numel(Time) > 1, sfreq = 1 / mean(diff(Time)); end

    ChannelFlag = i_get(rec, 'ChannelFlag', ones(nCh, 1));
    ChannelFlag = ChannelFlag(:);

    % Channel names/types travel WITH the store: the MEG-and-good selection downstream is
    % built from them, and having to reopen study.mat for it would undo the paging.
    ChannelName = {};  ChannelType = {};
    if isfield(S, 'chan') && isstruct(S.chan)
        ChannelName = i_get(S.chan, 'Name', {});
        ChannelType = i_get(S.chan, 'Type', {});
    end

    Comment = i_get(rec, 'Comment', '');
    Events  = i_get(rec, 'Events',  []);
    nAvg    = i_get(rec, 'nAvg',    1);
    source  = struct('studyDir', i_get(S, 'studyDir', ''), ...
                     'dataName', i_get(S, 'dataName', ''), ...
                     'study',    sfile);

    file = fullfile(dsdir, 'recording.mat');
    builtin('save', file, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', ...
        'ChannelType', 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', ...
        '-v7.3', '-nocompression');

    d = dir(file);
    fprintf('  %d ch x %d samples @ %g Hz (%.1f s) -> %s (%.0f MB)\n', ...
        nCh, nT, sfreq, nT / sfreq, file, d.bytes / 1e6);
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
