function fH = fbd_fig(doExport, pos)
% FBD_FIG  A demo figure, hidden when exporting so a headless run does not flash windows.
%   fH = fbd_fig(doExport, [x y w h])
% Author: Diellor Basha, 2026
    if nargin < 2 || isempty(pos), pos = [80 80 1100 700]; end
    % ⚠ AN EXPLICIT ROOT DEFAULT OF 'off' WINS. Setting a per-figure 'Visible' silently
    % overrides set(0,'DefaultFigureVisible','off'), so a harness that asks for headless
    % figures gets visible ones anyway. That is not academic: on a machine with no hardware
    % GL for MATLAB, a suite that runs these demos repeatedly creates hundreds of visible
    % figures and DEADLOCKS -- 0% CPU, no error, forever. Honour the default when it has
    % been set; otherwise behave exactly as before.
    vis = fbd_ternary(doExport, 'off', 'on');
    if strcmp(get(0, 'DefaultFigureVisible'), 'off'), vis = 'off'; end
    fH = figure('Color', 'w', 'Position', pos, 'Visible', vis);
end

% Author: Diellor Basha, 2026
