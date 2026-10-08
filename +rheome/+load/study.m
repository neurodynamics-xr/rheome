function st = study(name)
% LOAD.STUDY  Load the cached Brainstorm study structures from +data.
%
%   st = rheome.load.study(name)
%
% Returns a struct with the fields:
%   .chan  (rheome.io.read.channel)   .hm  (rheome.io.read.headmodel, leadfield)
%   .rec   (rheome.io.read.recording) .ncov (rheome.io.read.noisecov)
%   .studyDir, .dataName
% Written by rheome.import.study / rheome.import.dataset. Errors (with a hint) if not imported.
%
% See also: rheome.import.study, rheome.source.dirac, rheome.load.surface, rheome.load.dirac
%
% Author: Diellor Basha, 2026

    f = fullfile(rheome.load.root(), char(name), 'study.mat');
    if ~exist(f, 'file')
        error('load:study:missing', ...
            'No study cached for ''%s''. Run  rheome.import.study(''%s'', studyDir, dataName)  first.', name, name);
    end
    st = builtin('load', f);
end

% Author: Diellor Basha, 2026
