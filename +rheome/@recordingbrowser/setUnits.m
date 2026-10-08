function setUnits(app, ids)
% SETUNITS  Choose the channels stacked in the channels view.
%
%   app.setUnits(1:16)
%   app.setUnits(["MLO32" "MRO32"])
%
% The first one becomes the focused channel: it is what the raw load reads and what the
% detail view shows, so switching views never lands on a channel you did not pick.
%
% ⚠ MORE THAN MaxRows IS DRAWN AS MaxRows. Each row is a filled band of a few thousand
% vertices; a few dozen draw instantly and a few hundred do not. The extra rows stay in the
% selection and the readout says how many are shown.
%
% Author: Diellor Basha, 2026

    if isempty(ids)
        error('recordingbrowser:units', 'Select at least one channel.');
    end
    if isstring(ids) || iscellstr(ids) || ischar(ids)
        ids = cellstr(ids);
        u = arrayfun(@(i) rheome.recordingbrowser.unitid(app.Db, app.Scope, ids{i}), 1:numel(ids));
    else
        u = arrayfun(@(v) rheome.recordingbrowser.unitid(app.Db, app.Scope, v), double(ids(:)'));
    end
    app.Units = unique(u, 'stable');
    app.Unit  = app.Units(1);
    app.refresh();     % loaded samples stay with their channel, drawn on its own row
end
% Author: Diellor Basha, 2026
