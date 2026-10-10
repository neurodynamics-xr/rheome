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

    d = getenv(envvar);
    if ~isempty(d), ld_reachable(d, envvar); return; end
    here = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));   % toolbox root (parent of +rheome)
    if exist(fullfile(here, '.git'), 'file')                         % folder, or file in a worktree
        d = fullfile(here, inCheckout);
        ld_reachable(d, envvar);
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

function ld_reachable(d, envvar)
% A folder that cannot exist here is an error that says why, not a later 'file not found': the path is a
% link whose target is gone, or it lies on a drive (/Volumes/<name>) that is not mounted. A folder that
% merely does not exist yet (a fresh clone, a node-local cache) is fine -- the first import creates it.
    if isfolder(d), return; end
    vol = regexp(d, '^/Volumes/[^/]+', 'match', 'once');
    if ~isempty(vol) && ~isfolder(vol)
        error('load:noroot', ['rheome %s not found at %s: the drive %s is not mounted. Mount it, or set %s ' ...
            'to another folder (setenv(''%s'', ''/path'')).'], envvar, d, vol, envvar, envvar);
    end
    if usejava('jvm') && java.nio.file.Files.isSymbolicLink(java.io.File(d).toPath())
        error('load:noroot', ['rheome %s not found at %s: it is a link to %s, which is absent (is its drive ' ...
            'mounted?). Mount it, or set %s to another folder.'], envvar, d, ...
            char(java.nio.file.Files.readSymbolicLink(java.io.File(d).toPath())), envvar);
    end
end

% Author: Diellor Basha, 2026
