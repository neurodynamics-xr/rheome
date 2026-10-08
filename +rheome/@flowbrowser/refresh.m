function refresh(app)
% REFRESH  Redraw what depends on the current FRAME -- on the VISIBLE tab only.
%
% ⭐ Redrawing hidden tabs is pure waste and it is what makes stepping feel slow: the Cortex
% tab costs a GEMM over 20484 vertices per frame, the Spectrum tab a replot of 800 points.
% Paying for both when only one is on screen halves the frame time for nothing. Switching
% tabs redraws, so nothing is ever stale.
%
% The Timeline is drawn once per page -- only its cursor moves.
%
% Author: Diellor Basha, 2026

    tab = '';
    if ~isempty(app.Tabs) && isvalid(app.Tabs) && ~isempty(app.Tabs.SelectedTab)
        tab = app.Tabs.SelectedTab.Title;
    end
    switch tab
        case 'Cortex',   app.drawCortex();
        case 'Spectrum', app.drawSpectrum();
    end
    app.updateCursors();
end

% Author: Diellor Basha, 2026
