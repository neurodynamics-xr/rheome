function hAx = vectors(surf, V3, varargin)
% SHOW.VECTORS  Render a per-vertex 3-vector field on the cortex (magnitude + arrows).
%
%   rheome.show.vectors(surf, V3)
%   rheome.show.vectors(surf, V3, name, value, ...)
%
% Colors the surface by the vector magnitude and overlays quiver arrows at the
% strongest vertices, to visualize a reconstructed source current field.
%
% INPUTS:
%   surf : surface struct (rheome.io.read.surface) or file path
%   V3   : [nV x 3] vectors, or [3nV x 1] interleaved [x1,y1,z1, x2,...]
%
% OPTIONS (name/value):
%   'NumArrows' number of arrows to draw, at the top-magnitude vertices (default 400)
%   'Scale'     quiver length scale in mm-equivalent surface units (default auto)
%   'Color'     arrow color (default [0 0 0])
%   'View'      camera [az el] (default [-90 10])
%   'Parent'    axes to draw into (default new figure)
%   'Visible'   'on' (default) | 'off'
%   'Title'     title string
%
% OUTPUT:
%   hAx : axes handle.
%
% See also: rheome.show.surface, rheome.forward.reconstruct, rheome.inverse.dirac
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('NumArrows', 400);
    p.addParameter('Scale', []);
    p.addParameter('Color', [0.1 1 0.3]);
    p.addParameter('View', [-90 10]);
    p.addParameter('Parent', []);
    p.addParameter('Visible', 'on');
    p.addParameter('Title', []);
    p.addParameter('CLim', []);
    p.addParameter('Colorbar', true);
    p.addParameter('Colormap', hot(256));
    p.addParameter('Zoom', 1);         % camera zoom factor (>1 zooms in)
    p.addParameter('Target', []);      % [1x3] camera target (default: field centroid)
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(surf) || (isstring(surf) && isscalar(surf)), surf = rheome.io.read.surface(char(surf)); end
    V = surf.Vertices;  nV = size(V, 1);
    if isvector(V3) && numel(V3) == 3*nV
        V3 = [V3(1:3:end), V3(2:3:end), V3(3:3:end)];   % de-interleave
    end
    if size(V3,1) ~= nV || size(V3,2) ~= 3
        error('show:vectors:size', 'V3 must be [nV x 3] or [3nV x 1] (nV = %d).', nV);
    end

    mag = sqrt(sum(V3.^2, 2));

    % base surface colored by magnitude (sequential; magnitude is one-signed)
    hAx = rheome.show.surface(surf, mag, 'View', opt.View, 'Parent', opt.Parent, ...
                       'Visible', opt.Visible, 'Title', opt.Title, 'Colormap', opt.Colormap, ...
                       'CLim', opt.CLim, 'Colorbar', opt.Colorbar);

    % arrows at the strongest vertices
    nA = min(opt.NumArrows, nV);
    [~, ord] = sort(mag, 'descend');
    sel = ord(1:nA);
    if isempty(opt.Scale)
        bbox  = max(max(V, [], 1) - min(V, [], 1));
        scale = 0.06 * bbox / max(mag(sel) + eps);   % ~6% of the bounding box for the biggest arrow
    else
        scale = opt.Scale;
    end
    hold(hAx, 'on');
    quiver3(hAx, V(sel,1), V(sel,2), V(sel,3), ...
            scale*V3(sel,1), scale*V3(sel,2), scale*V3(sel,3), ...
            0, 'Color', opt.Color, 'LineWidth', 1.1, 'MaxHeadSize', 0.6);
    hold(hAx, 'off');

    % zoom the camera onto the active region so the arrows are readable
    if opt.Zoom ~= 1 || ~isempty(opt.Target)
        tgt = opt.Target;
        if isempty(tgt)
            w = mag / max(sum(mag), eps);  tgt = w' * V;    % magnitude-weighted field centroid
        end
        camtarget(hAx, tgt);
        if opt.Zoom ~= 1, camzoom(hAx, opt.Zoom); end
    end

    if nargout == 0, clear hAx; end
end

% Author: Diellor Basha, 2026
