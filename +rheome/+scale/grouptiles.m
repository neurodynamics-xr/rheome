function [P, A, ids, C] = grouptiles(S, Tp, h, l)
% SCALE.GROUPTILES  A participant's vertices in the group tiles of depth l, through the registered sphere.
%
%   [P, A, ids] = rheome.scale.grouptiles(S, Tp, "L", 3)
%   [P, A, ids, C] = rheome.scale.grouptiles(...)        % + tile centres
%
% Each participant vertex of hemisphere h belongs to the tile of its nearest template vertex on the sphere
% (Tp from rheome.scale.template, depth >= l), so tile k is the same anatomical place in every participant.
%   P    [nVh x nTiles] sparse 0/1 membership (rows = S.B.(h).gv order)
%   A    [nTiles x nTiles] logical adjacency (tiles sharing a mesh edge on the participant's cortex)
%   ids  [nTiles x 1] tile codes at depth l
%   C    struct: .u [nTiles x 3] tile centre on the template sphere (unit), .x [nTiles x 3] on the template
%        cortex (template coordinates, m; Brainstorm SCS: x anterior, y left, z up)
%
% See also: rheome.scale.template, rheome.scale.measure_atlas
%
% Author: Diellor Basha, 2026

    th = Tp.(h);  gv = S.B.(h).gv(:);  Fc = double(S.B.(h).S.Faces);
    us = double(S.Sf.Sphere(gv, :));  us = us ./ vecnorm(us, 2, 2);
    c = th.code(knnsearch(th.u, us));  sh = 2^(th.dmax - l);
    ids = unique(floor(th.code / sh));  [~, ti] = ismember(floor(c / sh), ids);
    P = sparse(1:numel(gv), ti, 1, numel(gv), numel(ids));
    E = [Fc(:,[1 2]); Fc(:,[2 3]); Fc(:,[3 1])];  x = ti(E(:,1)) ~= ti(E(:,2));
    A = sparse(ti(E(x,1)), ti(E(x,2)), 1, numel(ids), numel(ids));  A = (A + A') > 0;
    if nargout > 3
        [~, tt] = ismember(floor(th.code / sh), ids);  n = accumarray(tt, 1);
        u = [accumarray(tt, th.u(:,1)) accumarray(tt, th.u(:,2)) accumarray(tt, th.u(:,3))] ./ n;
        C = struct('u', u ./ vecnorm(u, 2, 2), ...
                   'x', [accumarray(tt, th.x(:,1)) accumarray(tt, th.x(:,2)) accumarray(tt, th.x(:,3))] ./ n);
    end
end

% Author: Diellor Basha, 2026
