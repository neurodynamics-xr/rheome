function X = harmonicfill(S, X, cortex)
% OPERATORS.HARMONICFILL  Replace a field on the medial wall by its harmonic extension from the cortical rim.
%
%   X = rheome.operators.harmonicfill(S, X, cortex)    % X [nV x m] or [3nV x m] (3 rows per vertex)
%   S.Res.ImagingKernel = rheome.operators.harmonicfill(S.Sf, S.Res.ImagingKernel, S.Sf.Cortex)
%
% The operators (LB eigenmodes, connection Laplacian, Dirac, div/curl, gauge) need the closed
% hemisphere: chi = 2, no boundary term. The wall is not cortex, so the field there is not data
% (rheome.geom.cortexmask). The fill solves Laplace's equation on the wall with the rim values as
% Dirichlet data -- the minimum-Dirichlet-energy completion -- so the operators see a closed field
% that adds no divergence or curl the rim does not force (review §4, layer 3).
%
%   L_ww u_w = -L_wc x_c        (L the cotan stiffness; w wall, c cortex)
%
% It is linear, so filling an inverse kernel's rows fills every source estimate made from it: no new
% inverse. A [3nV x m] input is an unconstrained kernel ordered (v-1)*3 + k; each orientation is
% filled on its own (the ambient field's x, y, z components).
%
% ⚠ The fill is a completion, not a measurement: anything read with its core on the wall (a
%   critical point, a vortex, a peak) is flagged and not reported (rheome.flow.critical .onWall).
% ⚠ Every wall component must touch cortex; a wall island with no rim has no Dirichlet data (error).
%
% See also: rheome.geom.cortexmask, rheome.operators.laplace_beltrami, rheome.scale.sensors
%
% Author: Diellor Basha, 2026

    nV = size(S.Vertices, 1);  cortex = logical(cortex(:));
    if numel(cortex) ~= nV, error('operators:harmonicfill:size', 'cortex has %d entries, the surface %d vertices', numel(cortex), nV); end
    w = find(~cortex);
    if isempty(w), return; end
    r = size(X, 1) / nV;
    if r ~= 1 && r ~= 3, error('operators:harmonicfill:rows', 'X needs nV or 3nV rows (%d vertices), got %d', nV, size(X, 1)); end
    L = rheome.operators.laplace_beltrami(double(S.Vertices), double(S.Faces));
    Lwc = L(w, cortex);
    [Rc, p] = chol(L(w, w));                              % SPD iff every wall component touches the rim
    if p, error('operators:harmonicfill:island', 'a medial-wall component has no cortical rim: no Dirichlet data'); end
    for k = 1:r
        rows = (0:nV-1) * r + k;
        X(rows(w), :) = -(Rc \ (Rc' \ (Lwc * X(rows(cortex), :))));
    end
end

% Author: Diellor Basha, 2026
