function key = meshkey(F)
% GEOM.MESHKEY  The editor's fingerprint of a mesh: FNV-1a over its faces' bytes, '<nF>:<hash>'.
%
%   key = rheome.geom.meshkey(F)        % F [nF x 3], 1-based
%
% The faces are hashed as the editor holds them: 0-based int32, row by row, little-endian. A stored product (a Hodge-Dirac
% basis, fused weights) records the key of the mesh it was made on; the editor refuses one whose key differs.
% Port of the Cortical Flow editor's meshKey (nxr-convert eigen/glue.py mesh_key).
%
% See also: rheome.geom.glue
%
% Author: Diellor Basha, 2026 (port of the Cortical Flow editor's meshKey)
bytes = typecast(int32(reshape((F - 1)', [], 1)), 'uint8');      % the hub is little-endian, as the editor's
h = uint64(2166136261);
for b = bytes'
    h = mod(bitxor(h, uint64(b)) * 16777619, 4294967296);          % < 2^56: exact in uint64
end
key = sprintf('%d:%08x', size(F,1), h);
end
