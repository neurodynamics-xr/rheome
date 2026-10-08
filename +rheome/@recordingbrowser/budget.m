function b = budget(app)
% BUDGET  The tile budget in force: an explicit number, or the axis width in pixels.
%
%   b = app.budget()
%
% ⭐ ONE TILE PER PIXEL IS THE RIGHT TARGET. Below it the view draws boxes where the store
% has detail; above it the extra tiles land on top of each other and cost work nobody can
% see. The axis knows its own width, so the default asks it rather than guessing, and a
% resized window changes the answer.
%
% ⚠ AN UNRENDERED FIGURE STILL HAS A POSITION. getpixelposition works with Visible off, so
% a headless script gets the figure's nominal width rather than a guess; the 1200 fallback
% is only for the case where there is no axis at all yet.
%
% Author: Diellor Basha, 2026

    b = app.Budget;
    if isnumeric(b) && isscalar(b) && isfinite(b) && b >= 1
        b = double(b);  return
    end
    b = 1200;
    if ~isempty(app.StatAx) && isvalid(app.StatAx)
        try
            p = getpixelposition(app.StatAx);
            if p(3) > 20, b = round(p(3)); end
        catch
        end
    end
end
% Author: Diellor Basha, 2026
