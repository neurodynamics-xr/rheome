function hFig = ridgestrip(surf, U, out, frames, varargin)
% SHOW.RIDGESTRIP  Filmstrip of a field with tracked ridges overlaid frame by frame.
%
%   rheome.show.ridgestrip(surf, U, out, frames)
%   hFig = rheome.show.ridgestrip(surf, U, out, frames, name, value, ...)
%
% The Lagrangian companion to rheome.show.filmstrip: each panel renders the field at one frame AND
% draws the motifs that are alive then -- a marker at the current position and a trail back to
% the ridge's birth. Reading left to right shows whether motifs persist in place or travel, on
% the data they were found in.
%
% INPUTS:
%   surf    surface struct (.Vertices, .Faces, .nV)
%   U       [nV x nT] the field the ridges were found in (|.| is taken)
%   out     output of rheome.detect.ridges
%   frames  vector of frame indices to render
%
% OPTIONS:
%   'Top'        show only the N longest-lived ridges (default 12; Inf = all)
%   'View'       camera [az el] (default [-90 10])
%   'CLim'       shared colour limits (default from U over the shown frames)
%   'Colormap'   background colormap (default parula)
%   'Lift'       marker/trail offset along the normal, fraction of mesh extent (default 0.012)
%   'TrailFrames' how many frames of trail to draw behind the current position (default Inf)
%   'Zoom'       camera zoom factor applied after rendering (default 1.8)
%   'Only'       restrict to these ridge indices before the Top sort
%   'Titles'     {1 x numel(frames)} panel titles (default 'frame k')
%   'Title'      overall title
%   'Visible'    'on' (default) | 'off'
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.detect.ridges, rheome.show.ridges, rheome.show.filmstrip
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Top', 12);
    p.addParameter('View', [-90 10]);
    p.addParameter('CLim', []);
    p.addParameter('Colormap', parula(256));
    p.addParameter('Lift', 0.012);
    p.addParameter('TrailFrames', Inf);
    p.addParameter('Zoom', 1.8);
    p.addParameter('Only', []);   % restrict to these ridge indices (before the Top sort)
    p.addParameter('Titles', {});
    p.addParameter('Title', 'tracked motifs');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    o = p.Results;

    R = out.ridges;
    if isempty(R), error('show:ridgestrip:empty', 'No ridges.'); end
    if ~isempty(o.Only), R = R(o.Only); end
    [~, ord] = sort([R.nFrames], 'descend');
    if isfinite(o.Top), ord = ord(1:min(o.Top, numel(ord))); end
    R = R(ord);

    U = abs(U);
    nF = numel(frames);
    if isempty(o.CLim), o.CLim = [min(U(:,frames),[],'all') max(U(:,frames),[],'all')]; end
    if isempty(o.Titles)
        o.Titles = arrayfun(@(k) sprintf('frame %d', k), frames, 'UniformOutput', false);
    end

    if isfield(surf,'VertNormals') && ~isempty(surf.VertNormals)
        Nrm = surf.VertNormals ./ max(vecnorm(surf.VertNormals,2,2), eps);
    else
        Nrm = surf.Vertices ./ max(vecnorm(surf.Vertices,2,2), eps);
    end
    span = max(max(surf.Vertices, [], 1) - min(surf.Vertices, [], 1));
    off  = o.Lift * span;

    cols = lines(max(numel(R), 7));
    hFig = figure('Color','w','Position',[40 40 300*nF 330],'Visible',o.Visible);

    for i = 1:nF
        f  = frames(i);
        ax = axes('Parent', hFig, 'Position', [(i-1)*0.955/nF + 0.004, 0.02, 0.955/nF - 0.006, 0.84]); %#ok<LAXES>
        rheome.show.surface(surf, U(:,f), 'Parent', ax, 'View', o.View, 'CLim', o.CLim, ...
            'Colormap', o.Colormap, 'Title', '');
        cb = findobj(hFig,'Type','colorbar');  delete(cb);      % one shared scale, drawn at the end
        camzoom(ax, o.Zoom);
        hold(ax,'on');
        for k = 1:numel(R)
            r   = R(k);
            hit = find(r.t == f, 1);
            if isempty(hit), continue; end
            lo  = max(1, hit - o.TrailFrames);
            P   = r.pos(lo:hit, :) + off * Nrm(r.v(lo:hit), :);
            plot3(ax, P(:,1), P(:,2), P(:,3), '-', 'LineWidth', 3, 'Color', cols(k,:));
            plot3(ax, P(end,1), P(end,2), P(end,3), 'o', 'MarkerSize', 15, ...
                'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', 'w', 'LineWidth', 1.6);
        end
        title(ax, o.Titles{i}, 'FontWeight','normal', 'FontSize', 10);
    end

    cbAx = axes('Parent', hFig, 'Position', [0.965 0.10 0.01 0.70], 'Visible','off');
    colormap(cbAx, o.Colormap);  caxis(cbAx, o.CLim);
    cb = colorbar(cbAx, 'Position', [0.968 0.14 0.008 0.62]);
    cb.Label.String = 'field amplitude';

    sgtitle(hFig, o.Title, 'FontWeight','bold', 'FontSize', 12);
end

% Author: Diellor Basha, 2026
