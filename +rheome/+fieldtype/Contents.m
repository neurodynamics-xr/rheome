% FIELDTYPE  What a field IS: the registry of field types, and the checks that use it.
%
% A field is a typed payload over a domain (+domain). Its type is the element rank it sits
% on, the bundle each element carries, the number field, and for a form its degree. Those
% decide the scalar count and decide which operators accept it.
%
%   rheome.fieldtype.registry    - every variant, as a table (ids mirror nxr-compute)
%   rheome.fieldtype.describe    - one variant by id, as a struct
%   rheome.fieldtype.components  - scalars per element (scalar 1, tangent 1, ambient 3, immersion 4)
%   rheome.fieldtype.validate    - does this array hold that type on that domain (scalar count)
%   rheome.fieldtype.matches     - are two types structurally the same
%
% ⚠ Spatial only: a [nV x nT] matrix is nT fields of one type, and time is its own domain.
%
% See also: rheome.domain.of, rheome.operators.registry, docs/OPERATOR-REGISTRY.md
%
% Author: Diellor Basha, 2026
