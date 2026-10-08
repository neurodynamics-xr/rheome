function hFig = filmstrip(surf, U, frames, varargin)
% SHOW.FILMSTRIP  Cortex snapshots of a space-time field at selected time frames.
%
%   rheome.show.filmstrip(surf, U, frames)
%   rheome.show.filmstrip(surf, U, frames, name, value, ...)
%
% Lays out a row of cortex renders (one per frame) with a SHARED colour scale, so the
% evolution of a space-time field (e.g. a joint-filter impulse response) can be read
% left-to-right. Can place its row into a multi-filter montage grid.
%
% INPUTS:
%   surf   surface struct or file path
%   U      [nV x nT] space-time field
%   frames vector of frame indices to show
%
% OPTIONS:
%   'Parent'   figure handle to draw into (default: new figure)
%   'Grid'     [nRows nCols] montage layout (default [1 numel(frames)])
%   'Row'      which row of the grid this filmstrip occupies (default 1)
%   'CLim'     shared colour limits (default symmetric from U(:,frames))
%   'Colormap' default: diverging (signed) via rheome.show.surface
%   'View'     camera [az el] (default [-90 10])
%   'Titles'   {1 x numel(frames)} per-panel titles (default '')
%   'RowLabel' text label for the row (shown on the first panel's y-axis)
%   'Visible'  'on' (default) | 'off'
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.filters.impulse, rheome.show.surface, rheome.show.gif
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Parent', []);
    p.addParameter('Grid', []);
    p.addParameter('Row', 1);
    p.addParameter('CLim', []);
    p.addParameter('Colormap', []);
    p.addParameter('View', [-90 10]);
    p.addParameter('Titles', []);
    p.addParameter('RowLabel', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(surf) || (isstring(surf) && isscalar(surf)), surf = rheome.io.read.surface(char(surf)); end
    nF = numel(frames);
    grid = opt.Grid;  if isempty(grid), grid = [1 nF]; end

    if isempty(opt.Parent)
        hFig = figure('Color','w','Visible',opt.Visible,'Position',[60 60 260*nF 260*grid(1)]);
    else
        hFig = opt.Parent;
    end

    clim = opt.CLim;
    if isempty(clim)
        % robust limits from the 99.5th percentile, so the initial delta spike does
        % not swamp the propagated field (which is the dynamics we want to see)
        vals = abs(reshape(U(:,frames), [], 1));
        q = quantile(vals, 0.995);  if q <= 0, q = max([vals; eps]); end
        if any(U(:,frames) < 0, 'all') && any(U(:,frames) > 0, 'all')
            clim = [-q q];                 % signed -> symmetric
        else
            clim = [max(min(U(:,frames),[],'all'), 0) q];
        end
    end

    for c = 1:nF
        ax = subplot(grid(1), grid(2), (opt.Row-1)*grid(2) + c, 'Parent', hFig);
        args = {'Parent',ax, 'View',opt.View, 'CLim',clim, 'Colorbar',false};
        if ~isempty(opt.Colormap), args = [args, {'Colormap',opt.Colormap}]; end %#ok<AGROW>
        if ~isempty(opt.Titles),   args = [args, {'Title',opt.Titles{c}}];     end %#ok<AGROW>
        rheome.show.surface(surf, U(:,frames(c)), args{:});
        if c == 1 && ~isempty(opt.RowLabel)
            ylabel(ax, opt.RowLabel, 'Visible','on', 'Rotation',90, ...
                   'FontWeight','bold', 'Interpreter','none');
        end
    end
    if nargout == 0, clear hFig; end
end

% Author: Diellor Basha, 2026
