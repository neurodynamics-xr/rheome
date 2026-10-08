function f = describe(id)
% FIELDTYPE.DESCRIBE  One row of the registry, as a struct, by id.
%
%   f = rheome.fieldtype.describe("ambientVertexWorld")
%
% Errors by name rather than returning empty: a typo in a field type should fail where it
% is written, not three operators later as a shape mismatch.
%
% Author: Diellor Basha, 2026

    T = rheome.fieldtype.registry();
    i = find(T.id == string(id), 1);
    if isempty(i)
        error('fieldtype:describe:unknown', ...
              'No field type ''%s''. rheome.fieldtype.registry() lists %d.', char(string(id)), height(T));
    end
    f = table2struct(T(i, :));
end
% Author: Diellor Basha, 2026
