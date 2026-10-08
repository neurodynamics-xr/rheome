function buildCortexAxes(app)
% BUILDCORTEXAXES  Create one surface per spatial scale, ONCE.
%
% Rebuilt only when the number of scales changes; drawCortex then only rewrites colour data.
%
% Author: Diellor Basha, 2026

    % panels: [1] all scales, then the dominant map IF enabled, then one per scale
    nS = app.Page.NumScales;
    nP = nS + 1 + app.ShowDominant;
    if ~isempty(app.CortexPatch) && numel(app.CortexPatch) == nP && all(isvalid(app.CortexPatch))
        return;
    end

    % ⚠ THE BANK'S MEMBER ORDER IS NOT SCALE ORDER: member 1 is the lowpass (127 mm) and the
    % wavelets then run FINE to COARSE (36 39 45 59 79 99 mm). Laying the panels out by
    % member index puts the scales on screen scrambled, so scanning left to right for "which
    % scale" reads a sequence that is not one. Panels are ordered COARSE -> FINE, matching
    % the scalogram's bottom -> top.
    kc = double(reshape(centerWavenumbers(app.Page.GraphBank), [], 1));
    [~, app.ScaleOrder] = sort(kc, 'ascend');            % small k = coarse, first
    % the hue ramp is assigned by DISPLAY position, so it is monotonic in actual scale --
    % otherwise the dominant-scale map's "warm coarse -> cool fine" would be neither.
    ramp = flipud(turbo(nS));
    app.ScaleColors = zeros(nS, 3);
    app.ScaleColors(app.ScaleOrder, :) = ramp;

    delete(app.CortexAx.Children);
    nCol = ceil(nP/2);
    app.CortexAx.RowHeight   = repmat({'1x'}, 1, 2);
    app.CortexAx.ColumnWidth = repmat({'1x'}, 1, nCol);

    S = app.Bundle.surface;
    app.CortexPatch = gobjects(1, nP);
    app.CortexTitle = gobjects(1, nP);
    for g = 0:nP-1
        pnl = uipanel(app.CortexAx, 'BorderType', 'none');
        gg  = uigridlayout(pnl, [2 1]);  gg.RowHeight = {18, '1x'};  gg.Padding = [0 0 0 0];
        app.CortexTitle(g+1) = uilabel(gg, 'Text', '', 'HorizontalAlignment', 'center', ...
            'Interpreter', 'tex', 'FontSize', 11);
        ax = uiaxes(gg);
        app.CortexPatch(g+1) = patch(ax, 'Faces', S.Faces, 'Vertices', S.Vertices, ...
            'FaceVertexCData', zeros(S.nV,1), 'FaceColor', 'interp', 'EdgeColor', 'none', ...
            'FaceLighting', 'gouraud', 'AmbientStrength', 0.55, 'DiffuseStrength', 0.5);
        axis(ax, 'equal', 'off', 'tight');
        view(ax, [-90 12]);                       % left lateral
        camlight(ax, 'headlight');
        ax.Interactions = [rotateInteraction zoomInteraction];
    end
end

% Author: Diellor Basha, 2026
