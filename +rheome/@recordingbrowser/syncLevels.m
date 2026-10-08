function syncLevels(app)
% SYNCLEVELS  Fill the level menu with the levels that actually exist, labelled by tile length.
%
% ⚠ A MENU ENTRY THAT DOES NOTHING IS A BUG. Below the diagonal's floor the store holds no
% rows, so picking such a level clamps straight back and the window does not move -- on the
% reference store in channel scope that was seven of thirteen entries, and the only sign was a
% line of grey text. The menu is rebuilt from the floor up instead, so everything in it
% changes something, and the entry says how long its tiles are.
%
% Rebuilt only when the floor or the level moves, since setting String resets the selection.
%
% Author: Diellor Basha, 2026

    if isempty(app.LevelSpin) || ~isvalid(app.LevelSpin), return; end
    g = app.Db.grid;
    ids = app.Floor:g.Lmax;
    if ~isequal(ids, app.LevelIds_)
        items = arrayfun(@(L) sprintf('L%d · %g s', L, g.tExtent(L+1)), ids, 'UniformOutput', false);
        app.LevelIds_ = ids;
        set(app.LevelSpin, 'Value', 1, 'String', items);   % Value first: a shorter list
    end                                                    % would leave it out of range
    i = find(app.LevelIds_ == app.Level, 1);
    if isempty(i), i = 1; end
    set(app.LevelSpin, 'Value', i);
end
% Author: Diellor Basha, 2026
