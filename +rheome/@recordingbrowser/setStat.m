function setStat(app, name)
% SETSTAT  Choose the statistic drawn on the top axis, from rheome.select.derive's vocabulary.
%
% A per-band statistic (bandPower, share) is drawn for the selected Band and raises the
% floor to that band's natural level; a time or spectral statistic does not.
%
% Author: Diellor Basha, 2026

    name = char(name);
    V = rheome.select.derive();
    if ~ismember(string(name), V.name)
        error('recordingbrowser:stat', 'Unknown statistic ''%s''. rheome.select.derive() lists the vocabulary.', name);
    end
    ok = rheome.recordingbrowser.statsfor(app.Db);
    if ~ismember(string(name), ok)
        error('recordingbrowser:stat', ...
            'This store cannot serve ''%s''; it offers %s.', name, strjoin(cellstr(ok), ', '));
    end
    app.Stat = name;
    app.refresh();
end
% Author: Diellor Basha, 2026
