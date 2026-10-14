function finish(F, ttl, stem)
% FINISH  Close a house figure: its title goes to the figure name (a journal figure's title is the first
% sentence of its legend, never in the artwork), F is kept for later rows and nxr_save, and with a file
% stem the figure is exported at print size (PDF with fonts embedded, SVG, PNG at 600 dpi).
%
% Author: Diellor Basha, 2026
    if ~isempty(ttl), F.fig.Name = char(ttl); end
    setappdata(F.fig, 'nxr', F);
    if ~isempty(stem), nxr_save(F, char(stem)); end
end
