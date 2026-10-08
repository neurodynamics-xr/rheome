function L = levelfor(db, width, budget, floorLevel)
% RECORDINGBROWSER.LEVELFOR  The level-of-detail rule: the finest level under the tile budget.
%
%   L = rheome.recordingbrowser.levelfor(db, width, budget, floorLevel)
%
% A window `width` seconds wide holds ceil(width / tExtent(L)) tiles at level L, which halves
% with every coarser level. Take the FINEST level whose count is at or under `budget`, then
% clamp to [floorLevel, Lmax]: that keeps the cost of a redraw bounded no matter how far out
% the view is, and it is the only thing that decides which level a zoom reads.
%
% ⚠ THE FLOOR WINS. Below floorLevel the store carries nothing for the unit or the band
% (the diagonal), so zooming further cannot buy resolution -- it narrows the window until
% the raw axis takes over instead.
%
% Author: Diellor Basha, 2026

    g = db.grid;
    L = max(0, min(floorLevel, g.Lmax));
    while L < g.Lmax && ceil(width / g.tExtent(L+1)) > budget
        L = L + 1;
    end
end
% Author: Diellor Basha, 2026
