function v = nearest_vertices(templateSphere, templateFaces, templateHemi, pointsSphere, pointsHemi)
% REGISTER.NEAREST_VERTICES  Points (cores, sources, peaks) given on the registered sphere, each to the template vertex it falls on.
%
%   v = rheome.register.nearest_vertices(templateSphere, templateFaces, templateHemi, pointsSphere, pointsHemi)     % v [n x 1], 1-based
%
% The template's sphere sampled AT the point, nearest (FreeSurfer's vertex correspondence). An error names how many points
% reach no template vertex. As nxr.register.vertices.nearest_vertices.
%
% See also: rheome.register.registration, rheome.register.register_values
%
% Author: Diellor Basha, 2026
arguments, templateSphere (:,3) double, templateFaces (:,3) double, templateHemi (:,1) char, pointsSphere (:,3) double, pointsHemi (:,1) char, end
R = rheome.register.registration(templateSphere, templateFaces, templateHemi, pointsSphere, pointsHemi, Sampling="nearest");
c = full(sum(R ~= 0, 2));
if any(c ~= 1), error('rheome:register:points', '%d points reach no template vertex', nnz(c ~= 1)); end
[~, v] = max(R, [], 2); v = full(v);
end
