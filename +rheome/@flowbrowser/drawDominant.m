function drawDominant(app, Z)
% DRAWDOMINANT  One cortex coloured by WHICH spatial scale dominates at each vertex.
%
% ⭐ WHAT SEPARATE PANELS CANNOT SHOW. Per-scale surfaces separate exactly what you want to
% compare. To ask "is the flow coarse here and fine there" the cortex must be held fixed
% while scale varies WITHIN it. Hue = the winning scale (warm coarse, cool fine); brightness
% = how much energy is there at all.
%
% OFF BY DEFAULT (rheome.flowbrowser.ShowDominant). Nine panels crowd the Cortex tab; this is kept
% live rather than deleted because it answers a question the per-scale panels structurally
% cannot -- set ShowDominant = true to bring it back.
%
% ⚠ BRIGHTNESS IS NOT DECORATION. An argmax over a flat spectrum is meaningless and would
% flicker frame to frame; darkening low-energy vertices toward the background is what stops
% the map asserting a scale where there is nothing to assert. The gamma is a display choice
% and errs toward showing less.
%
% Author: Diellor Basha, 2026

    fp = app.Page;
    if nargin < 2, Z = []; end
    [iDom, total] = dominantScale(fp, app.Frame, Z);

    hue = app.ScaleColors(iDom, :);                       % [V x 3]
    w   = total / max(max(total), realmin);
    w   = w .^ 0.45;                                      % gamma: lift mid energies
    bg  = [0.72 0.72 0.74];                               % neutral cortex, not black
    rgb = bg + w .* (hue - bg);

    if ~app.ShowDominant, return; end
    h = app.CortexPatch(2);
    if isvalid(h), set(h, 'FaceVertexCData', rgb); end
end

% Author: Diellor Basha, 2026
