function name = rheomeTestSubject(kind)
% RHEOMETESTSUBJECT  The cached subject the data-dependent tests read; unset, those tests are skipped.
%
%   name = rheomeTestSubject()           % $RHEOME_TEST_SUBJECT: a subject imported with rheome.import.*
%   name = rheomeTestSubject('second')   % $RHEOME_TEST_SUBJECT2: a second anatomy
%   name = rheomeTestSubject('noise')    % $RHEOME_TEST_NOISE: an empty-room study for the same sensors
%
% Most tests build their own synthetic data and need none of this. The few that check a result on a
% real cortex read the cache under rheome.load.root(); without the variable they name a subject that
% is not there, fail to load it and are filtered by assumption, never failed.
%
% See also: rheome.load.root, rheome.import.surface, rheome.import.study
%
% Author: Diellor Basha, 2026

    if nargin < 1, kind = 'subject'; end
    var = struct('subject', 'RHEOME_TEST_SUBJECT', 'second', 'RHEOME_TEST_SUBJECT2', 'noise', 'RHEOME_TEST_NOISE');
    name = getenv(var.(kind));
    if isempty(name), name = ['rheome_test_' kind '_unset']; end
end

% Author: Diellor Basha, 2026
