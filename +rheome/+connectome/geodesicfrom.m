function d = geodesicfrom(S, member, source)
% CONNECTOME.GEODESICFROM  Geodesic distance on the cortical surface from a source parcel's centroid, per parcel.
%
%   d = rheome.connectome.geodesicfrom(S, member, source)
%
% For H5/TP5 (D_region): the source point is the source parcel's vertex nearest its Euclidean
% centroid; the distance (heat method, rheome.geom.geodesic -- the port of nxr-compute's
% geodesic.cpp) is computed on the connected surface component holding that vertex, and each
% parcel gets the mean over its vertices. Parcels on another component (the other hemisphere)
% are NaN: a geodesic does not cross the midline.
%
% INPUTS:
%   S       surface with .Vertices [nV x 3] and .Faces [nF x 3] (e.g. the subject's white or mid surface)
%   member  [nP x nV] logical/sparse parcel membership (rheome.io.read.atlas .Membership)
%   source  index of the source parcel (e.g. entorhinal L)
%
% OUTPUT:
%   d [nP x 1] mean geodesic distance from the source centroid (mesh units), NaN off-component
%
% Author: Diellor Basha, 2026

    V = double(S.Vertices);  F = double(S.Faces);  nV = size(V, 1);
    srcVerts = find(member(source, :));
    [~, i] = min(sum((V(srcVerts, :) - mean(V(srcVerts, :), 1)).^2, 2));
    seed = srcVerts(i);

    A = sparse(F, F(:, [2 3 1]), 1, nV, nV);
    comp = conncomp(graph((A + A') > 0));
    onComp = (comp == comp(seed))';
    keepFace = all(onComp(F), 2);
    newIndex = zeros(nV, 1);  newIndex(onComp) = 1:nnz(onComp);
    sub = struct('Vertices', V(onComp, :), 'Faces', newIndex(F(keepFace, :)));

    dv = zeros(nV, 1);
    dv(onComp) = rheome.geom.geodesic(sub, newIndex(seed));
    W = double(member);  on = double(onComp);
    d = (W * dv) ./ (W * on);                           % mean over on-component vertices
    d(W * on == 0) = NaN;                               % parcel wholly on another component
end

% Author: Diellor Basha, 2026
