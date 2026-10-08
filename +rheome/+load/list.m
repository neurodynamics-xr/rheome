function names = list()
% LOAD.LIST  Print the cached datasets in +data and what each holds.
%
%   rheome.load.list
%   names = rheome.load.list()
%
% See also: rheome.import.dataset, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    d = rheome.load.root();
    if ~exist(d, 'dir')
        fprintf('No +data cache yet (%s).\n', d);
        if nargout > 0, names = {}; end
        return;
    end
    e = dir(d);  e = e([e.isdir] & ~startsWith({e.name}, '.'));
    fprintf('Cached datasets in %s:\n', d);
    out = {};
    for i = 1:numel(e)
        nm = e(i).name;
        f  = dir(fullfile(d, nm, '*.mat'));
        parts = {};
        if any(strcmp({f.name}, 'surface.mat')), parts{end+1} = 'surface'; end
        if any(strcmp({f.name}, 'study.mat')),   parts{end+1} = 'study';   end
        dl = f(startsWith({f.name}, 'dirac__'));
        for j = 1:numel(dl), parts{end+1} = strrep(strrep(dl(j).name, '.mat',''), 'dirac__', 'dirac '); end %#ok<AGROW>
        fprintf('  %-16s [%s]  (%.0f MB)\n', nm, strjoin(parts, ', '), sum([f.bytes])/1e6);
        out{end+1,1} = nm; %#ok<AGROW>
    end
    if isempty(e), fprintf('  (none)\n'); end
    if nargout > 0, names = out; end
end

% Author: Diellor Basha, 2026
