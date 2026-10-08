function file = study(name, studyDir, dataName)
% IMPORT.STUDY  Read the Brainstorm study structures + cache them to +data.
%
%   file = rheome.import.study(name, studyDir, dataName)
%
% Reads the channel file, the unconstrained MEG leadfield (headmodel_surf_os_meg*),
% the recording, and the noise covariance (rheome.io.read.*), and writes them to
% +data/<name>/study.mat. This caches the (large) leadfield/recording locally so the
% source demos never re-read the external volume.
%
% dataName may be a bare filename (joined to studyDir) OR an ABSOLUTE PATH. The absolute form
% is what you need when the recording and the geometry live in different Brainstorm conditions
% -- e.g. a continuous recording imported into its own condition while the channel file and
% leadfield stay in the original @raw condition, which is the normal layout for a raw study
% that was imported afterwards.
%
% See also: rheome.import.surface, rheome.import.dirac, rheome.import.dataset, rheome.load.study, rheome.source.dirac
%
% Author: Diellor Basha, 2026

    dsdir = fullfile(rheome.load.root(), char(name));
    if ~exist(dsdir, 'dir'), mkdir(dsdir); end

    chanF = i_find(studyDir, 'channel_*.mat');
    hmF   = i_find(studyDir, 'headmodel_surf_os_meg*.mat');   % newest match, ignores ._ sidecars
    ncovF = fullfile(studyDir, 'noisecov_full.mat');
    dataName = char(dataName);
    if startsWith(dataName, filesep) || ~isempty(regexp(dataName, '^[A-Za-z]:', 'once'))
        dataF = dataName;                       % absolute: recording lives in another condition
    else
        dataF = fullfile(studyDir, dataName);
    end
    for f = {chanF, hmF, ncovF, dataF}
        if isempty(f{1}) || ~exist(f{1}, 'file')
            error('import:study:missing', 'Missing study file: %s', f{1});
        end
    end

    fprintf('rheome.import.study[%s]: reading channel / leadfield / recording / noisecov\n', name);
    chan = rheome.io.read.channel(chanF);       %#ok<NASGU>
    hm   = rheome.io.read.headmodel(hmF);
    rec  = rheome.io.read.recording(dataF);     %#ok<NASGU>
    ncov = rheome.io.read.noisecov(ncovF);      %#ok<NASGU>

    file = fullfile(dsdir, 'study.mat');
    builtin('save', file, 'chan', 'hm', 'rec', 'ncov', 'studyDir', 'dataName', '-v7.3');
    fprintf('  cached study (leadfield %dx%d) -> %s\n', size(hm.Gain,1), size(hm.Gain,2), file);
end

% newest matching file, ignoring macOS ._ sidecars (mirrors rheome.source.dirac's i_find)
function f = i_find(d, pat)
    L = dir(fullfile(d, pat));  f = '';
    L = L(~startsWith({L.name}, '._'));
    if isempty(L), return; end
    [~, newest] = max([L.datenum]);
    f = fullfile(d, L(newest).name);
end

% Author: Diellor Basha, 2026
