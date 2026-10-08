function fbd_save(fH, outDir, name)
% FBD_SAVE  Export a figure as a PNG when an output folder was given.
% Author: Diellor Basha, 2026
    if isempty(outDir), return; end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    exportgraphics(fH, fullfile(outDir, [name '.png']), 'Resolution', 130);
end

% Author: Diellor Basha, 2026
