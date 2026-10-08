function s = statsfor(db)
% RECORDINGBROWSER.STATSFOR  The statistics a store can actually serve, as names.
%
% A preview store holds min and max only, so offering rms in its menu would be offering an
% error. Everything else serves the whole vocabulary.
%
% Author: Diellor Basha, 2026

    V = rheome.select.derive();
    if isfield(db, 'preview') && db.preview
        s = intersect(V.name, ["min","max","ptp"], 'stable');
    else
        s = V.name(:)';
    end
end
% Author: Diellor Basha, 2026
