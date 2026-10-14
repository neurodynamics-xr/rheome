function finish(F, ttl, stem)
% FINISH  Close a house figure: its title goes to the figure name (a journal figure's title is the first
% sentence of its legend, never in the artwork), every text holding $...$ is set by LaTeX (tokens type.math:
% symbols and unit exponents are maths, never Unicode glyphs), the descriptive names of the colormaps used
% go to getappdata(F.fig, 'colormapNames') for the legend (tokens colormapName, never 'PuOr_r'), F is kept
% for later rows and nxr_save, and with a file stem the figure is exported at print size.
%
% Author: Diellor Basha, 2026
    if ~isempty(ttl), F.fig.Name = char(ttl); end
    h = findall(F.fig, '-property', 'Interpreter');
    for i = 1:numel(h)
        if isprop(h(i), 'String'), s = h(i).String; elseif isprop(h(i), 'Label'), s = h(i).Label; else, continue; end
        if iscell(s), s = strjoin(s, ' '); end
        if contains(string(s), '$'), h(i).Interpreter = 'latex'; end
    end
    % shortcut: nxr_colorbar sets its key-band label in the DEFAULT theme's face; put plain text in this theme's
    % face until the kit takes the theme from F (asked of Graphics)
    t = findall(F.fig, 'Type', 'text', 'Interpreter', 'none');  set(t, 'FontName', F.tokens.fontName);
    cb = findall(F.fig, 'Type', 'ColorBar');  q = {};
    for i = 1:numel(cb)
        p = strsplit(cb(i).Tag, '|');  if numel(p) >= 2, q{end+1} = p{2}; end %#ok<AGROW>
    end
    q = unique(q, 'stable');
    setappdata(F.fig, 'colormapNames', cellfun(@nxr_cmap_name, q, 'UniformOutput', false));
    setappdata(F.fig, 'nxr', F);
    if ~isempty(stem), nxr_save(F, char(stem)); end
end
