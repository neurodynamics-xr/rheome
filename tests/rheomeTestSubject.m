function name = rheomeTestSubject(varargin)
% RHEOMETESTSUBJECT  The cached subject the data-dependent tests read; unset, those tests are skipped.
%
%   name = rheomeTestSubject()           % $RHEOME_TEST_SUBJECT: a subject imported with rheome.import.*
%   name = rheomeTestSubject('second')   % $RHEOME_TEST_SUBJECT2: a second anatomy
%   name = rheomeTestSubject('noise')    % $RHEOME_TEST_NOISE: an empty-room study for the same sensors
%   name = rheomeTestSubject(tc, kinds...)  % the same, and ASSUME each kind is set and cached, saying why not
%
% Most tests build their own synthetic data and need none of this. The few that check a result on a
% real cortex read the cache under rheome.load.root(). Call the tc form once, before the test's own
% try/catch (a bare catch would swallow the assumption and its message): the test is then filtered
% with the reason -- variable unset, cache folder absent (drive not mounted), or subject not cached --
% never failed.
%
% See also: rheome.load.root, rheome.import.surface, rheome.import.study
%
% Author: Diellor Basha, 2026

    tc = [];
    if nargin && isa(varargin{1}, 'matlab.unittest.TestCase'), tc = varargin{1}; varargin(1) = []; end
    kinds = string(varargin);
    if isempty(kinds), kinds = "subject"; end
    var = struct('subject', 'RHEOME_TEST_SUBJECT', 'second', 'RHEOME_TEST_SUBJECT2', 'noise', 'RHEOME_TEST_NOISE');
    what = struct('subject', 'subject', 'second', 'second subject', 'noise', 'empty-room study');
    names = strings(size(kinds));
    for k = 1:numel(kinds)
        names(k) = getenv(var.(kinds(k)));
        if ~isempty(tc)
            tc.assumeNotEmpty(char(names(k)), sprintf('%s is unset: no cached %s to test on', var.(kinds(k)), what.(kinds(k))));
            try
                root = rheome.load.root();
            catch err
                tc.assumeFail(err.message);           % the drive or the link is gone: load.root says which
            end
            tc.assumeTrue(isfolder(root), sprintf('rheome +data cache not found at %s: set RHEOME_DATA', root));
            tc.assumeTrue(isfolder(fullfile(root, names(k))), ...
                sprintf('%s = %s is not in the +data cache at %s', var.(kinds(k)), names(k), root));
        end
        if names(k) == "", names(k) = "rheome_test_" + kinds(k) + "_unset"; end
    end
    name = char(names(1));
end

% Author: Diellor Basha, 2026
