function h = hemispheres(P, F, left)
% GEOM.HEMISPHERES  'L' / 'R' for each vertex of a two-piece cortex: its two largest connected pieces, told apart by an axis.
%
%   h = rheome.geom.hemispheres(P, F, "+y")      % P [n x 3], F [nF x 3] 1-based; h [n x 1] char, '?' off the two pieces
%
% left  the axis on which the left hemisphere lies: "+y" (Brainstorm SCS), "-x" (MNI, fsaverage); the piece whose mean position
%       along it is larger is 'L'.
% Port of nxr-convert eigen/register.py hemispheres.
%
% See also: rheome.geom.sphereregister, rheome.eigen.transfer
%
% Author: Diellor Basha, 2026 (port of nxr-convert hemispheres)
arguments, P (:,3) double, F (:,3) double, left (1,1) string {mustBeMember(left, ["+x" "-x" "+y" "-y" "+z" "-z"])}, end
n = size(P, 1);
sgn = 1 - 2 * startsWith(left, "-"); ax = find('xyz' == char(extractAfter(left, 1)));
comp = conncomp(graph(sparse(F(:), reshape(F(:, [2 3 1]), [], 1), 1, n, n) + sparse(reshape(F(:, [2 3 1]), [], 1), F(:), 1, n, n)))';
[~, big] = sort(accumarray(comp, 1), 'descend'); big = big(1:2);
side = [mean(P(comp == big(1), ax)) mean(P(comp == big(2), ax))] * sgn;
h = repmat('?', n, 1);
[~, l] = max(side); h(comp == big(l)) = 'L'; h(comp == big(3 - l)) = 'R';
end
