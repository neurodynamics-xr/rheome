function tf = has(name)
% LOAD.HAS  True if 'name' is a cached dataset (a folder under +data), not a file path.
%
%   tf = rheome.load.has(name)
%
% Used by the demos / rheome.source.dirac to tell a dataset NAME ('mysubject') from a file PATH
% ('/…/tess_cortex.mat'): a dataset name has no path separators, no .mat extension, and
% resolves to a +data/<name>/ folder.
%
% See also: rheome.load.dataset, rheome.import.dataset
%
% Author: Diellor Basha, 2026

    tf = ischar(name) && ~isempty(name) ...
         && ~any(name == '/' | name == filesep) ...
         && ~endsWith(lower(name), '.mat') ...
         && exist(fullfile(rheome.load.root(), name), 'dir') == 7;
end

% Author: Diellor Basha, 2026
