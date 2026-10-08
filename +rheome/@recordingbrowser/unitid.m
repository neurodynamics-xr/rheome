function u = unitid(db, scope, want)
% RECORDINGBROWSER.UNITID  A channel name, a channel id or a group id -> the unit id for a scope.
%
%   u = rheome.recordingbrowser.unitid(db, 'channel', 'MLO32')
%   u = rheome.recordingbrowser.unitid(db, 'group', [])          % the root of the sensor tree
%
% Empty picks a default: the first good channel, or the tree root, so the browser always
% opens on something that exists.
%
% Author: Diellor Basha, 2026

    if strcmp(scope, 'group')
        nodes = db.groupNodes;
        if isempty(nodes)
            error('recordingbrowser:groups', ...
                'This store has no sensor groups (built without positions); Scope must be ''channel''.');
        end
        if isempty(want), u = min(nodes); return; end        % node 1 is the root
        if ischar(want) || isstring(want)
            error('recordingbrowser:unit', 'Group scope takes a node id, not a name.');
        end
        u = double(want);
        if ~ismember(u, nodes)
            error('recordingbrowser:unit', 'Node %d carries no rows. Internal nodes are %s.', u, mat2str(nodes(:)'));
        end
        return
    end
    names = string(db.meta.ChannelName(:));
    if isempty(want), u = 1; return; end
    if ischar(want) || isstring(want)
        i = find(names == string(want), 1);
        if isempty(i)
            error('recordingbrowser:unit', 'No channel named ''%s'' in this store.', char(string(want)));
        end
        u = i;  return
    end
    u = double(want);
    if u < 1 || u > numel(names) || mod(u, 1) ~= 0
        error('recordingbrowser:unit', 'Channel id must be an integer in 1..%d.', numel(names));
    end
end
% Author: Diellor Basha, 2026
