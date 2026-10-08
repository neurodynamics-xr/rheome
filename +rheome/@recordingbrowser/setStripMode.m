function setStripMode(app, m)
% SETSTRIPMODE  'natural' (each band at its own tile length) or 'level' (all at one level).
%
% 'natural' shows the constant-Q staircase and fills the whole frequency range, because
% every band is drawn at a level that carries it. 'level' lines the strip's columns up with
% the statistic above, and shades the bands the level does not hold.
%
% Author: Diellor Basha, 2026

    m = char(m);
    if ~ismember(m, {'natural','level'})
        error('recordingbrowser:stripmode', 'StripMode must be ''natural'' or ''level'', got ''%s''.', m);
    end
    app.StripMode = m;
    app.refresh();
end
% Author: Diellor Basha, 2026
