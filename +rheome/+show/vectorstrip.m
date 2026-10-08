function hFig = vectorstrip(surf, J, frames, varargin)
% SHOW.VECTORSTRIP  Cortex quiver snapshots of a space-time VECTOR field at time frames.
%
%   rheome.show.vectorstrip(surf, J, frames)
%   rheome.show.vectorstrip(surf, J, frames, name, value, ...)
%
% The vector-field analogue of rheome.show.filmstrip: a row of cortex renders (one per frame)
% showing the current DIRECTION as arrows and the magnitude as colour, with a SHARED
% arrow scale and colour range so the evolution of an unconstrained (3-vector) source
% reads left-to-right. Can place its row into a multi-filter montage grid.
%
% INPUTS:
%   surf   surface struct or file path
%   J      [3nV x nT] vector field, rows [x1,y1,z1, x2,y2,z2, ...]
%   frames vector of frame indices to show
%
% OPTIONS:
%   'Parent'    figure handle to draw into (default: new figure)
%   'Grid'      [nRows nCols] montage layout (default [1 numel(frames)])
%   'Row'       which grid row this strip occupies (default 1)
%   'NumArrows' arrows per panel, at the top-magnitude vertices (default 200)
%   'View'      camera [az el] (default [-90 10])
%   'Colormap'  magnitude colormap (default hot)
%   'RowLabel'  text label for the row
%   'Titles'    {1 x numel(frames)} per-panel titles
%   'Visible'   'on' (default) | 'off'
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.show.vectors, rheome.show.filmstrip, rheome.filters.impulse, rheome.forward.reconstruct
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Parent', []);
    p.addParameter('Grid', []);
    p.addParameter('Row', 1);
    p.addParameter('NumArrows', 200);
    p.addParameter('View', [-90 10]);
    p.addParameter('Colormap', hot(256));
    p.addParameter('RowLabel', '');
    p.addParameter('Titles', []);
    p.addParameter('Visible', 'on');
    p.addParameter('Zoom', 2.4);
    p.addParameter('Target', []);
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(surf) || (isstring(surf) && isscalar(surf)), surf = rheome.io.read.surface(char(surf)); end
    V  = surf.Vertices;  nV = size(V,1);
    nF = numel(frames);
    grid = opt.Grid;  if isempty(grid), grid = [1 nF]; end

    if isempty(opt.Parent)
        hFig = figure('Color','w','Visible',opt.Visible,'Position',[60 60 300*nF 300*grid(1)]);
    else
        hFig = opt.Parent;
    end

    % shared magnitude scale + arrow scale from the strongest vertex across the shown frames
    Jf   = J(:, frames);
    magF = sqrt( Jf(1:3:end,:).^2 + Jf(2:3:end,:).^2 + Jf(3:3:end,:).^2 );   % [nV x nF]
    gmax = max(magF(:));  if gmax <= 0, gmax = eps; end
    arrowScale = 0.07 * max(max(V, [], 1) - min(V, [], 1)) / gmax;    % biggest arrow (any frame) ~7% of the bbox

    % shared camera target: magnitude-weighted centroid over all shown frames (so the
    % zoomed view is stable across the strip)
    tgt = opt.Target;
    if isempty(tgt)
        w = sum(magF,2);  w = w / max(sum(w), eps);  tgt = w' * V;
    end

    for c = 1:nF
        ax = subplot(grid(1), grid(2), (opt.Row-1)*grid(2) + c, 'Parent', hFig);
        V3 = [J(1:3:end,frames(c)), J(2:3:end,frames(c)), J(3:3:end,frames(c))];
        ttl = '';  if ~isempty(opt.Titles), ttl = opt.Titles{c}; end
        rheome.show.vectors(surf, V3, 'Parent',ax, 'View',opt.View, 'NumArrows',opt.NumArrows, ...
            'Scale',arrowScale, 'CLim',[0 gmax], 'Colorbar',false, 'Colormap',opt.Colormap, ...
            'Title',ttl, 'Zoom',opt.Zoom, 'Target',tgt);
        if c == 1 && ~isempty(opt.RowLabel)
            ylabel(ax, opt.RowLabel, 'Visible','on', 'Rotation',90, 'FontWeight','bold', 'Interpreter','none');
        end
    end
    if nargout == 0, clear hFig; end
end

% Author: Diellor Basha, 2026
