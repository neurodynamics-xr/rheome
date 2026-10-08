function hAx = surface(surf, data, varargin)
% SHOW.SURFACE  Render a cortical surface, optionally colored by a per-vertex field.
%
%   rheome.show.surface(surf)                 % plain gray cortex
%   rheome.show.surface(surf, data)           % cortex colored by a [nV x 1] scalar field
%   rheome.show.surface(file, ...)            % 'file' = path to a Brainstorm surface .mat
%   hAx = rheome.show.surface(..., name,value)
%
% The workhorse viewer for the demo: it draws the mesh once and, when given a
% per-vertex scalar (an eigenmode, a heat kernel, a realized atom's magnitude),
% shades the surface by it. Reused by every later visualization.
%
% INPUTS:
%   surf : struct from rheome.io.read.surface, OR a path to a Brainstorm surface .mat,
%          OR a [nV x 3] vertex matrix (then pass faces via the 'Faces' option).
%   data : optional [nV x 1] per-vertex scalar to color the surface. [] = flat gray.
%
% OPTIONS (name/value):
%   'Faces'      [nF x 3] faces, required only if 'surf' is a bare vertex matrix
%   'Colormap'   colormap for the overlay (default: diverging for signed data
%                (blue-white-red), sequential (parula) for one-signed data)
%   'FaceAlpha'  surface opacity 0..1 (default 1)
%   'EdgeColor'  edge color (default 'none')
%   'View'       [az el] camera view (default [-90 10], left-lateral-ish)
%   'CLim'       [lo hi] color limits for the overlay (default auto, symmetric if signed)
%   'Parent'     axes handle to draw into (default: new figure + axes)
%   'Visible'    'on' (default) | 'off'  (use 'off' for headless rendering/export)
%   'Title'      title string (default: the surface Comment)
%
% OUTPUT:
%   hAx : handle to the axes the surface was drawn into.
%
% See also: rheome.io.read.surface, rheome.operators.laplace_beltrami
%
% Author: Diellor Basha, 2026

    if nargin < 2, data = []; end

    % --- resolve the surface argument into V, F ---
    p = inputParser;
    p.addParameter('Faces', []);
    p.addParameter('Colormap', []);
    p.addParameter('FaceAlpha', 1);
    p.addParameter('EdgeColor', 'none');
    p.addParameter('View', [-90 10]);
    p.addParameter('CLim', []);
    p.addParameter('Parent', []);
    p.addParameter('Visible', 'on');
    p.addParameter('Title', []);
    p.addParameter('Colorbar', true);
    p.parse(varargin{:});
    opt = p.Results;

    titleUnset = any(strcmp('Title', p.UsingDefaults));   % distinguish unset from explicit ''
    ttl = opt.Title;
    if ischar(surf) || (isstring(surf) && isscalar(surf))
        surf = rheome.io.read.surface(char(surf));
    end
    if isstruct(surf)
        V = surf.Vertices;  F = surf.Faces;
        if titleUnset && isfield(surf, 'Comment'), ttl = surf.Comment; end
    else
        V = surf;  F = opt.Faces;
        if isempty(F)
            error('show:surface:faces', 'When ''surf'' is a vertex matrix you must pass ''Faces''.');
        end
    end
    if ~isempty(data)
        data = data(:);
        if numel(data) ~= size(V, 1)
            error('show:surface:data', 'data has %d entries but the surface has %d vertices.', ...
                numel(data), size(V, 1));
        end
    end

    % --- axes ---
    if isempty(opt.Parent)
        hFig = figure('Color', 'w', 'Visible', opt.Visible);
        hAx  = axes('Parent', hFig);
    else
        hAx  = opt.Parent;
    end

    % --- draw ---
    if isempty(data)
        patch('Parent', hAx, 'Faces', F, 'Vertices', V, ...
              'FaceColor', [0.82 0.82 0.82], 'EdgeColor', opt.EdgeColor, ...
              'FaceAlpha', opt.FaceAlpha, 'FaceLighting', 'gouraud');
    else
        patch('Parent', hAx, 'Faces', F, 'Vertices', V, ...
              'FaceVertexCData', data, 'FaceColor', 'interp', ...
              'EdgeColor', opt.EdgeColor, 'FaceAlpha', opt.FaceAlpha, ...
              'FaceLighting', 'gouraud');
        isSigned = any(data < 0) && any(data > 0);
        cmap = opt.Colormap;
        if isempty(cmap)
            if isSigned, cmap = i_diverging(256); else, cmap = parula(256); end
        end
        colormap(hAx, cmap);
        if isempty(opt.CLim)
            if isSigned                            % signed field -> symmetric limits
                m = max(abs(data));  if m > 0, set(hAx, 'CLim', [-m m]); end
            else
                lo = min(data); hi = max(data);
                if hi > lo, set(hAx, 'CLim', [lo hi]); end
            end
        else
            set(hAx, 'CLim', opt.CLim);
        end
        if opt.Colorbar, colorbar(hAx); end
    end

    % --- presentation ---
    axis(hAx, 'equal');  axis(hAx, 'off');
    view(hAx, opt.View);
    camlight(hAx, 'headlight');  material(hAx, 'dull');
    if ~isempty(ttl) && (ischar(ttl) || isstring(ttl))
        title(hAx, ttl, 'Interpreter', 'none');
        set(get(hAx, 'Title'), 'Visible', 'on');   % keep title despite axis off
    end

    if nargout == 0, clear hAx; end
end

function cmap = i_diverging(n)
% Blue-white-red diverging colormap for signed fields (0 -> white).
    if nargin < 1, n = 256; end
    x = linspace(0, 1, n)';
    lo = [0.23 0.30 0.75];   % blue
    mid = [1 1 1];            % white
    hi = [0.75 0.15 0.15];   % red
    cmap = zeros(n, 3);
    for c = 1:3
        cmap(:,c) = interp1([0 0.5 1], [lo(c) mid(c) hi(c)], x, 'linear');
    end
end

% Author: Diellor Basha, 2026
