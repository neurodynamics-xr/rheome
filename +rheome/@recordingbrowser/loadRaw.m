function n = loadRaw(app)
% LOADRAW  Read the raw samples of the current window, on request.
%
%   n = app.loadRaw()      -> samples read (0 if it refused)
%
% ⭐ DELIBERATELY IMPERATIVE. Everything else in the browser is answered from the pyramid,
% which is why a zoom costs one array read; raw samples are the one thing that touches the
% recording, and on this store a read pulls EVERY channel's hyperslab, so 8 s of one sensor
% moves 10 MB and 2 minutes moves 156 MB. Reading that whenever a window happened to get
% short would make panning unpredictable, so the samples arrive when you ask and stay until
% you ask again.
%
% ⚠ IT REFUSES A WIDE WINDOW rather than stalling, and the panel says the size it would
% have moved. Raise RawSeconds if you mean it.
%
% What is loaded stays loaded, and it stays WITH ITS CHANNEL: the samples carry the unit
% they came from, so browsing to another row draws them on their own row rather than
% relabelling them, and the panel says whose they are.
%
% Author: Diellor Basha, 2026

    n = 0;
    if strcmp(app.Scope, 'group')
        app.setStatus('raw samples are per channel: switch Scope to channel');
        return
    end
    meta = app.Db.meta;
    w = app.Window;
    if diff(w) > app.RawSeconds
        app.setStatus(sprintf('refused: %.0f s window is wider than RawSeconds (%g s), %.0f MB off disk', ...
            diff(w), app.RawSeconds, diff(w) * meta.fs * meta.C * 8 / 1e6));
        return
    end
    if ~isfield(meta, 'source') || exist(char(meta.source), 'file') ~= 2
        app.setStatus('the recording store this was built from is not on this disk');
        return
    end
    try
        row = meta.iChannel(app.Unit);
        if isempty(app.Pr_) || app.PrUnit_ ~= row
            app.Pr_ = rheome.pagedrecording(char(meta.source), 'Channels', row, 'Precision', 'double');
            app.PrUnit_ = row;
        end
        a = max(1, floor(w(1) * meta.fs) + 1);
        b = min(meta.nT, ceil(w(2) * meta.fs));
        t = tic;
        x = read(app.Pr_, a, b).';
        app.Raw = struct('x', x, 't', ((a:b)' - 1) / meta.fs, 'span', [(a-1)/meta.fs, b/meta.fs], ...
                         'unit', app.Unit, 'note', '');
        n = numel(x);
        app.refresh();
        app.setStatus(sprintf('loaded %d samples (%.1f-%.1f s) in %.0f ms', n, (a-1)/meta.fs, b/meta.fs, 1e3*toc(t)));
    catch err
        app.setStatus(sprintf('raw read failed: %s', err.message));
    end
end
% Author: Diellor Basha, 2026
