function s = describe(id)
% SELECTION.DESCRIBE  One selection by id, as a struct.
%
%   s = rheome.selection.describe("cortex_tile")
%
% See also: rheome.selection.registry, rheome.selection.check
%
% Author: Diellor Basha, 2026
    T = rheome.selection.registry();
    j = T.id == string(id);
    if ~any(j)
        error('selection:describe:unknown', 'No selection "%s". rheome.selection.registry lists them.', id);
    end
    s = table2struct(T(j, :));
end

% Author: Diellor Basha, 2026
