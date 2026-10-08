function gif(surf, U, file, varargin)
% SHOW.GIF  Write an animated GIF of a space-time field evolving on the cortex.
%
%   rheome.show.gif(surf, U, file)
%   rheome.show.gif(surf, U, file, name, value, ...)
%
% Renders each time frame of U on the cortex (fixed colour scale and camera) and writes
% an animated GIF -- the true motion of a joint-filter impulse response or a filtered
% dynamics segment.
%
% INPUTS:
%   surf  surface struct or file path
%   U     [nV x nT] space-time field
%   file  output .gif path
%
% OPTIONS:
%   'Frames'    frame indices (default 1:nT)
%   'CLim'      colour limits (default symmetric from U(:,Frames))
%   'Colormap'  default: diverging (signed) / parula (one-signed)
%   'View'      camera [az el] (default [-90 10])
%   'DelayTime' seconds per frame (default 0.06)
%   'Resolution' export DPI (default 80)
%
% See also: rheome.show.filmstrip, rheome.filters.impulse
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Frames', []);
    p.addParameter('CLim', []);
    p.addParameter('Colormap', []);
    p.addParameter('View', [-90 10]);
    p.addParameter('DelayTime', 0.06);
    p.addParameter('Resolution', 80);
    p.parse(varargin{:});
    opt = p.Results;

    if ischar(surf) || (isstring(surf) && isscalar(surf)), surf = rheome.io.read.surface(char(surf)); end
    frames = opt.Frames;  if isempty(frames), frames = 1:size(U,2); end

    clim = opt.CLim;
    if isempty(clim)
        m = max(abs(reshape(U(:,frames), [], 1)));
        if any(U(:,frames) < 0, 'all') && any(U(:,frames) > 0, 'all'), clim = [-m m];
        else, clim = [min(U(:,frames),[],'all') max(U(:,frames),[],'all')]; end
        if clim(1) == clim(2), clim = clim + [-1 1]*eps; end
    end

    % Draw once; update the patch CData per frame and capture with getframe (fast).
    hFig = figure('Color','w','Visible','off','Position',[80 80 420 380]);
    ax   = axes('Parent', hFig);
    hAx  = rheome.show.surface(surf, U(:,frames(1)), 'Parent',ax, 'View',opt.View, 'CLim',clim, 'Colorbar',false, 'Title','');
    if ~isempty(opt.Colormap), colormap(hAx, opt.Colormap); end
    hp = findobj(hAx, 'Type','patch');

    for k = 1:numel(frames)
        set(hp, 'FaceVertexCData', U(:,frames(k)));
        drawnow;
        fr = getframe(hFig);
        [A, map] = rgb2ind(fr.cdata, 256);
        if k == 1
            imwrite(A, map, file, 'gif', 'LoopCount', Inf, 'DelayTime', opt.DelayTime);
        else
            imwrite(A, map, file, 'gif', 'WriteMode', 'append', 'DelayTime', opt.DelayTime);
        end
    end
    close(hFig);
end

% Author: Diellor Basha, 2026
