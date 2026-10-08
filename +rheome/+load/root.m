function d = root()
% LOAD.ROOT  Absolute path to the +data cache folder (shared by rheome.import.* and rheome.load.*).
%
%   d = rheome.load.root()          % $RHEOME_DATA if set, else <repo>/+data in a git clone,
%                            % else <userpath>/rheome/data for the installed toolbox
%
% Set RHEOME_DATA (setenv, or the shell) to cache a subject somewhere else -- on a cluster node
% that is node-local disk, so a per-subject job imports into node-local disk and never writes into
% the checkout. Mirrors rheome.load.outroot / RHEOME_OUT for the output side.
%
% See also: rheome.load.list, rheome.load.outroot, rheome.import.dataset, rheome.scale.run
%
% Author: Diellor Basha, 2026

    d = ld_home('RHEOME_DATA', '+data', 'data');
end

% Author: Diellor Basha, 2026
