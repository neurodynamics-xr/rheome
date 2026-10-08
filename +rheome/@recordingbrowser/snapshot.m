function f = snapshot(app, file)
% SNAPSHOT  Save the current view as a PNG -- the figure a paper or a note keeps.
%
%   f = app.snapshot()                       % _figures/recordingbrowser_<name>_L<level>.png
%   f = app.snapshot('scratch/view.png')
%
% An interactive view that cannot be exported leaves nothing behind; this is what makes the
% browser usable as a figure source and not only as a toy.
%
% ⚠ DRAWNOW FIRST, ALWAYS. An invisible figure has not rendered, and exportgraphics on one
% fails with "Figure must contain graphics" -- which is exactly the case a headless script or
% a test hits. Flushing the queue is the whole fix.
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(file)
        d = '_figures';
        if exist(d, 'dir') ~= 7, mkdir(d); end
        file = fullfile(d, sprintf('recordingbrowser_%s_%s_L%02d.png', app.Name, app.Stat, app.Level));
    else
        d = fileparts(file);
        if ~isempty(d) && exist(d, 'dir') ~= 7, mkdir(d); end
    end
    if isempty(app.Fig) || ~isvalid(app.Fig)
        error('recordingbrowser:snapshot', 'The window is closed.');
    end
    drawnow;
    exportgraphics(app.Fig, file, 'Resolution', 110);
    f = file;
    app.setStatus(sprintf('saved %s', file));
end
% Author: Diellor Basha, 2026
