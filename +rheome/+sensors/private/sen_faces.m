function F = sen_faces(arr)
% SEN_FACES  The readout triangulation -- faces (2-cells), and nothing else.
%
%   F = sen_faces(arr)
%
% ⭐ THIS IS NOT THE OPERATOR. Filtering runs on the Gaussian-weighted kNN graph and never
% touches these faces. They exist only to supply the faces (2-cells) that a winding number
% and a per-face wavevector require -- rheome.detect.phasesingularity and rheome.flow.phasegradient.
% Keeping the two apart is why the operator can stay purely GSPBox.
%
% ⭐ Dim 1 RETURNS EMPTY, AND THAT IS THE POINT. A winding number is summed around a face (a
% 2-cell), and a chain admits none, so a winding number is undefined on a probe rather than
% merely noisy. The emptiness IS the statement; nothing downstream needs a special case for
% probes.
%
% ⚠ convhull OF A CAP CLOSES IT. Projecting a helmet or scalp cap onto its fitted sphere and
% taking the hull yields a CLOSED surface, and the faces that close it span the opening --
% huge triangles joining sensors on opposite rims. Left in, they assert loops the array does
% not have. They are removed by edge length: any face with an edge beyond 3x the median.
%
% ⚠ THE AXIS-DROP (planar branch) IS NOT AN ORTHONORMAL PROJECTION IN GENERAL. Rotating into
% the SVD basis -- Q = Pc*V(:,1:2) -- was tried first and rejected: on an axis-aligned
% lattice, V's entries come back ULP-close to 0/+-1 but not bit-exact, and multiplying
% introduces ~1e-19-scale noise. A perfect grid's unit cells are exactly cocircular, so that
% noise is enough to flip qhull's in-circle tie-break per cell (and its boundary-point
% tolerance), inflating a 5x6 grid's face count from the analytic 40 to 50 (18 hull points
% become 7). Dropping the coordinate axis most aligned with the normal instead keeps Q
% bit-identical to the raw input for every planar array this project actually constructs
% (rheome.sensors.grid emits Pos(:,3) == 0 exactly). For a plane genuinely tilted relative to all
% three axes, the drop is an oblique projection, not a metric-preserving one: its worst case,
% when the normal is equally inclined to all three axes, scales areas by 1/sqrt(3) (~0.577),
% so triangle SHAPES can differ from a true in-plane Delaunay triangulation. It leaves V, E
% and F -- hence the topology, i.e. the faces (2-cells) -- unchanged, which is the only thing
% a winding number or a per-face wavevector consumes. Do not revert this to Pc*V(:,1:2); that
% reintroduces the qhull degeneracy above.
%
% INPUTS:
%   arr  a rheome.sensors.* array struct
% OUTPUT:
%   F    [nF x 3] vertex indices, or [] when arr.Dim == 1
%
% See also: rheome.sensors.graph, rheome.detect.phasesingularity, rheome.flow.phasegradient
%
% Author: Diellor Basha, 2026

    if arr.Dim < 2
        F = [];  return;
    end

    P = arr.Pos;
    Pc = P - mean(P, 1);
    [~, S, V] = svd(Pc, 'econ');
    sv = diag(S);
    isPlanar = sv(3) <= 1e-9 * sv(1);

    if isPlanar
        % Drop the coordinate axis most aligned with the plane normal -- see the AXIS-DROP
        % header note above for why this replaces a rotation into the SVD basis.
        [~, dropAxis] = max(abs(V(:, 3)));
        Q = Pc(:, setdiff(1:3, dropAxis));
        F = delaunay(Q(:,1), Q(:,2));
    else
        [c, R] = i_spherefit(P);
        U  = P - c.';
        rn = vecnorm(U, 2, 2);
        if any(rn <= 1e-9 * R)
            error('sensors:faces:sensorAtSphereCentre', ...
                ['A sensor lies within 1e-9*R (R = %.3g m) of the fitted sphere centre. ' ...
                 'Its radial direction is numerically undefined, and convhull would ' ...
                 'silently accept whatever direction floating-point noise happens to ' ...
                 'produce -- so this is refused rather than triangulated.'], R);
        end
        U = U ./ rn;                              % project onto the fitted sphere
        try
            F = convhull(U(:,1), U(:,2), U(:,3));
        catch me
            error('sensors:faces:hullFailed', ...
                'Could not triangulate the readout hull for %d sensors: %s', ...
                size(U, 1), me.message);
        end
        F = i_dropspan(P, F, 3.0);
    end
end

function [c, R] = i_spherefit(P)
% Least-squares sphere centre: |p|^2 = 2 p.c + (R^2 - |c|^2) is LINEAR in (c, R^2-|c|^2).
    A = [2*P, ones(size(P,1), 1)];
    b = sum(P.^2, 2);
    x = A \ b;
    c = x(1:3);
    R = sqrt(x(4) + sum(c.^2));
end

function F = i_dropspan(P, F, factor)
% Drop faces with any edge beyond factor x the median edge -- the cap-closing triangles.
    e = [vecnorm(P(F(:,1),:) - P(F(:,2),:), 2, 2), ...
         vecnorm(P(F(:,2),:) - P(F(:,3),:), 2, 2), ...
         vecnorm(P(F(:,3),:) - P(F(:,1),:), 2, 2)];
    F = F(max(e, [], 2) <= factor * median(e(:)), :);
end

% Author: Diellor Basha, 2026
