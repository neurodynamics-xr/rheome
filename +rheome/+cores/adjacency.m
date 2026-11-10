function A = adjacency(F, nv)
% CORES.ADJACENCY  Symmetric vertex adjacency of a triangle mesh (nxr support_adjacency, whole mesh as support).
%
%   A = rheome.cores.adjacency(F, nv)   % F [nF x 3] 1-based faces; A [nv x nv] sparse logical
%
% Ports nxr-analysis src/nxr/vortex.py support_adjacency (commit 903fdc2): vertices i,j adjacent when they share a face edge.
%
% See also: rheome.cores.patches
%
% Author: Diellor Basha, 2026
    F = double(F);
    A = sparse(F(:), reshape(circshift(F, 1, 2), [], 1), true, nv, nv);
    A = (A | A.');
end
