function setView(app, v)
% SETVIEW  'channels' (the stack of min/max envelopes) or 'detail' (one channel in full).
%
% ⚠ THE DETAIL VIEW NEEDS A TILE STORE. Its band strip and its statistics do not exist in a
% preview store, so asking for it there is refused by name rather than drawing empty panes.
%
% Author: Diellor Basha, 2026

    v = char(v);
    if ~ismember(v, {'channels','detail'})
        error('recordingbrowser:view', 'View must be ''channels'' or ''detail'', got ''%s''.', v);
    end
    if strcmp(v, 'detail') && isempty(app.Db.bands) && ~app.Envelope
        app.Envelope = true;            % a preview store can only serve the envelope
    end
    app.View = v;
    app.layout();
    app.refresh();
end
% Author: Diellor Basha, 2026
