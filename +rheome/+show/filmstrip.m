function hFig = filmstrip(surf, U, frames, varargin)
% SHOW.FILMSTRIP  Cortex snapshots of a space-time field at selected time frames.
%
%   rheome.show.filmstrip(surf, U, frames)
%   rheome.show.filmstrip(surf, U, frames, name, value, ...)
%
% Lays out a row of cortex renders (one per frame) with a SHARED colour scale, so the
% evolution of a space-time field (e.g. a joint-filter impulse response) can be read
% left-to-right. Can place its row into a multi-filter montage grid. The figure is journal width
% (190 mm), one panel per row with one colourbar under it; sizes and colours are the figure tokens.
%
% INPUTS:
%   surf   surface struct or file path
%   U      [nV x nT] space-time field
%   frames vector of frame indices to show
%
% OPTIONS:
%   'Parent'   a figure made by a strip before, to draw another row into (default: new figure)
%   'Grid'     [nRows nCols] montage layout (default [1 numel(frames)])
%   'Row'      which row of the grid this filmstrip occupies (default 1)
%   'CLim'     shared colour limits (default symmetric from U(:,frames))
%   'Colormap' default: the house map, flow-signed-field (signed) or flow-magnitude
%   'Quantity', 'Unit'  the colourbar label, 'quantity (unit)' (default 'field', 'a.u.')
%   'View'     camera [az el] (default [-90 10])
%   'Titles'   {1 x numel(frames)} per-panel titles (default '')
%   'RowLabel' the row's panel title, after its letter
%   'Title'    figure title: the figure's name (the journal legend carries it, not the artwork)
%   'Export'   file stem: write <stem>.pdf/.svg/.png at print size (nxr_save)
%   'Visible'  'on' (default) | 'off'
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.filters.impulse, rheome.show.surface, rheome.show.gif
%
% Author: Diellor Basha, 2026

    tokens();
    p = inputParser;
    p.addParameter('Parent', []);
    p.addParameter('Grid', []);
    p.addParameter('Row', 1);
    p.addParameter('CLim', []);
    p.addParameter('Colormap', []);
    p.addParameter('View', [-90 10]);
    p.addParameter('Titles', []);
    p.addParameter('RowLabel', '');
    p.addParameter('Quantity', 'field');
    p.addParameter('Unit', 'a.u.');
    p.addParameter('Title', '');
    p.addParameter('Export', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(surf) || (isstring(surf) && isscalar(surf)), surf = rheome.io.read.surface(char(surf)); end
    nF = numel(frames);
    grid = opt.Grid;  if isempty(grid), grid = [1 nF]; end

    clim = opt.CLim;
    signed = any(U(:,frames) < 0, 'all') && any(U(:,frames) > 0, 'all');
    if isempty(clim)
        % robust limits from the 99.5th percentile, so the initial delta spike does
        % not swamp the propagated field (which is the dynamics we want to see)
        vals = abs(reshape(U(:,frames), [], 1));
        q = quantile(vals, 0.995);  if q <= 0, q = max([vals; eps]); end
        if signed
            clim = [-q q];                 % signed -> symmetric
        else
            clim = [max(min(U(:,frames),[],'all'), 0) q];
        end
    end
    cmap = opt.Colormap;
    if isempty(cmap), cmap = nxr_cmap(ternary(signed, 'flow-signed-field', 'flow-magnitude')); end
    kind = ternary(signed && abs(clim(1) + clim(2)) <= 1e-9*max(abs(clim)), 'signed', 'magnitude');

    [F, axs, clim] = striprow(opt.Parent, grid, opt.Row, nF, opt.Visible, opt.RowLabel, ...
                        {opt.Quantity, opt.Unit, kind, clim, cmap});
    for c = 1:nF
        args = {'Parent',axs(c), 'View',opt.View, 'CLim',clim, 'Colorbar',false, 'Colormap',cmap, 'Title',''};
        if ~isempty(opt.Titles), args{end} = opt.Titles{c}; end
        rheome.show.surface(surf, U(:,frames(c)), args{:});
    end
    hFig = F.fig;
    finish(F, opt.Title, opt.Export);
    if nargout == 0, clear hFig; end
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
