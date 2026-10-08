function n = components(id)
% FIELDTYPE.COMPONENTS  Scalars per element of a field type.
%
%   n = rheome.fieldtype.components("immersionVertex")   -> 4
%
% scalar 1, tangent 1 (one complex coordinate per element), ambient 3, immersion 4. The
% only arithmetic a shape check needs, kept in one place so it cannot drift between the
% registry and the validator.
%
% Author: Diellor Basha, 2026

    if isstruct(id), n = id.components; return; end
    n = rheome.fieldtype.describe(id).components;
end
% Author: Diellor Basha, 2026
