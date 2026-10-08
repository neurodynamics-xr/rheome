function t = measure(app, kind, value, opts)
% MEASURE  Write a number onto the selected tile, keyed by the pyramid.
%
%   app.selectTime(265);  app.measure("amplitude", 4.2e-13)
%   app.measure("vortex", 1, Band=6, Source="manual")
%
% ⭐ THE KEY IS THE PYRAMID'S OWN. A measurement is addressed by (recording, scope, unit,
% level, k) and optionally a band -- the same key the moments use -- so the number lands on
% the tile you clicked, for the channel you were looking at, and joins straight back to the
% statistics of that tile and of any ancestor (rheome.select.context).
%
% ⚠ IT DOES NOT MERGE, AND IT IS NOT MEANT TO. This goes into the second bookkeeping matrix
% (rheome.select.measure, the `measure` relation, kind 'note'), which is sparse and per node: a
% value here says nothing about the parent tile, no query prunes on it, and the same kind at
% two levels is two observations. That is the freedom that makes it useful for anything a
% moment cannot express.
%
% Writes are transactions: atomic, idempotent by content hash, removed by rheome.select.rollback.
%
% See also: rheome.select.measure, rheome.select.measures, rheome.select.context, rheome.recordingbrowser/measurements
%
% Author: Diellor Basha, 2026

    arguments
        app
        kind  (1,1) string
        value (1,1) double
        opts.Band   (1,1) double = 0
        opts.Source (1,1) string = "browser"
        opts.Author (1,1) string = ""
    end
    s = app.Selected;
    if isempty(s) || ~isfield(s, 'k')
        error('recordingbrowser:measure', ...
            'No tile is selected. Click the view, or call selectTime(t), first.');
    end
    rows = table(s.level, s.k, value, s.unit, opts.Band, ...
                 'VariableNames', {'level','k','value','unit_id','band_id'});
    t = rheome.select.measure(app.Db, rows, Kind=kind, Scope=string(s.scope), ...
                       Source=opts.Source, Author=opts.Author);
    app.refreshMarks();
    app.setStatus(sprintf('%s = %.4g on %s tile (L%d, k%d), txn %d%s', kind, value, ...
        i_unit(app, s), s.level, s.k, t.txn_id, i_existed(t)));
end

function u = i_unit(app, s)
    if strcmp(s.scope, 'group'), u = sprintf('node %d', s.unit);
    else, u = char(string(app.Db.meta.ChannelName{s.unit})); end
end

function e = i_existed(t)
    if t.existed, e = ' (already written)'; else, e = ''; end
end
% Author: Diellor Basha, 2026
