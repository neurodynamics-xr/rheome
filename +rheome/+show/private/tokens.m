function T = tokens(theme)
% TOKENS  The lab's figure tokens (library/styles/tokens/figure-tokens.json) as a struct, and the house
% figure kit (nxr_figure, nxr_axes_style, nxr_colorbar, nxr_cmap, nxr_save) on the path.
%   The styles folder is $NXR_STYLES, else ~/workspace/library/styles. Every size and colour of rheome.show
%   comes from here: the figure rules are tokens read by every exporter, not literals (Diellor, 2026-10-14).
%   theme  'pergammon' (default) | 'neuroimage' (MS1) | 'nature': page widths, height and type family.
%
% Author: Diellor Basha, 2026
    persistent ready
    if isempty(ready)
        d = getenv('NXR_STYLES');
        if isempty(d), d = fullfile(getenv('HOME'), 'workspace', 'library', 'styles'); end
        m = fullfile(d, 'matlab');
        if ~isfile(fullfile(m, 'nxr_figure_tokens.m'))
            error('rheome:show:tokens', 'figure tokens not found in %s (set NXR_STYLES to library/styles)', m);
        end
        if ~contains([pathsep path pathsep], [pathsep m pathsep]), addpath(m); end
        ready = true;
    end
    if nargin < 1, theme = ''; end
    T = nxr_figure_tokens(theme);
end
