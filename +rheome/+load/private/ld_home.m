function d = ld_home(envvar, inCheckout, inUser)
% LD_HOME  Default folder for rheome.load.root / rheome.load.outroot: the env var, else beside the code in a git
% checkout, else under userpath for an installed add-on.
%
%   d = ld_home('RHEOME_DATA', '+data', 'data')
%
% ⚠ An installed .mltbx lives in the add-ons folder, which MATLAB may replace on upgrade or
% uninstall, so nothing is ever cached or written there. A clone (it has .git) keeps the old layout.
%
% See also: rheome.load.root, rheome.load.outroot
%
% Author: Diellor Basha, 2026

    d = getenv(envvar);
    if ~isempty(d), return; end
    here = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % toolbox root (parent of +rheome)
    if exist(fullfile(here, '.git'), 'file')                         % folder, or file in a worktree
        d = fullfile(here, inCheckout);
        return
    end
    u = userpath;
    if isempty(u)
        error('load:nohome', ['userpath is not set, so there is no default folder for %s. ' ...
            'Set it, e.g. setenv(''%s'', ''/path/to/folder''), or restore userpath (userpath(''reset'')).'], ...
            inUser, envvar);
    end
    d = fullfile(u, 'rheome', inUser);
end

% Author: Diellor Basha, 2026
