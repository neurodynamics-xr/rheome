function setScope(app, s)
% SETSCOPE  Switch between single channels and sensor-tree nodes (the position pyramid).
%
% Switching resets the unit to a default that exists in the new scope: a node id is not a
% channel id, and carrying one across would read a different sensor than the label says.
%
% Author: Diellor Basha, 2026

    s = char(s);
    if ~ismember(s, {'channel','group'})
        error('recordingbrowser:scope', 'Scope must be ''channel'' or ''group'', got ''%s''.', s);
    end
    app.Scope = s;
    app.Unit  = rheome.recordingbrowser.unitid(app.Db, s, []);
    app.Raw   = struct('t', [], 'x', [], 'span', [NaN NaN], 'note', '');
    app.Ladder = table();                      % the floor columns are scope-dependent
    app.syncUnits();
    app.refresh();
end
% Author: Diellor Basha, 2026
