function tf = matches(a, b)
% FIELDTYPE.MATCHES  Are two field types structurally the same?
%
%   tf = rheome.fieldtype.matches("scalarVertex", f)
%
% Structural means rank, bundle, number field, form degree and representation. The gauge
% and the n-RoSy order are advisory in the C++ registry and are not match keys here either,
% for the same reason: they parameterise an operator, they do not change what the payload is.
%
% Author: Diellor Basha, 2026

    if ~isstruct(a), a = rheome.fieldtype.describe(a); end
    if ~isstruct(b), b = rheome.fieldtype.describe(b); end
    tf = a.rank == b.rank && a.bundle == b.bundle && a.field == b.field && ...
         a.n_form == b.n_form && a.repr == b.repr;
end
% Author: Diellor Basha, 2026
