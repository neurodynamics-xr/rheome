function syncUnits(app)
% SYNCUNITS  Fill the unit menu for the current scope (channel names, or sensor-tree nodes).
%
% A popupmenu carries only strings, so the ids live beside it in UnitIds_ and the callback
% maps the selected index back to an id. Rebuilt on a scope change, not per redraw: a
% 270-channel list is cheap once and wasteful every frame.
%
% Author: Diellor Basha, 2026

    if strcmp(app.Scope, 'group')
        nodes = app.Db.groupNodes(:)';
        tr = app.Db.tree;
        items = cell(numel(nodes), 1);
        for i = 1:numel(nodes)
            r = find(tr.node_id == nodes(i), 1);
            items{i} = sprintf('node %d (depth %d, %d sensors)', nodes(i), tr.depth(r), tr.n_sensors(r));
        end
        app.UnitIds_ = nodes;
    else
        items = cellstr(string(app.Db.meta.ChannelName(:)));
        app.UnitIds_ = 1:numel(items);
    end
    app.UnitNames_ = items;
    if ~isempty(app.UnitDrop) && isvalid(app.UnitDrop)
        i = find(app.UnitIds_ == app.Unit, 1);
        if isempty(i), i = 1; end
        set(app.UnitDrop, 'String', items, 'Value', i);
    end
    if ~isempty(app.ChanList) && isvalid(app.ChanList)
        sel = find(ismember(app.UnitIds_, app.Units));
        if isempty(sel), sel = 1; end
        set(app.ChanList, 'Value', 1, 'String', items);
        set(app.ChanList, 'Value', sel);
    end
end
% Author: Diellor Basha, 2026
