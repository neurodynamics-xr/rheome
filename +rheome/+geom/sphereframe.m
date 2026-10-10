function fr = sphereframe(V, F, N, Sph, opts)
% GEOM.SPHEREFRAME  The group gauge: "north" and "west" on a cortex, from its FreeSurfer sphere.
%
%   fr = rheome.geom.sphereframe(V, F, N, Sph)                % poles at the registered sphere's +-z
%   fr = rheome.geom.sphereframe(V, F, N, Sph, Axis=[1 0 0])  % the second chart, poles at +-x
%
% THE GAUGE FOR GROUP ANALYSIS (Diellor, 2026-10-08). Its two singularities sit at the poles of
% FreeSurfer's registered sphere (sphere.reg, Brainstorm's Reg.Sphere), i.e. at the SAME anatomical
% place in every subject (fsaverage +z: precentral; -z: fusiform / inferior temporal), and e1 is the
% sphere's meridian towards +Axis PUSHED FORWARD to the cortex by the sphere->cortex vertex map (a
% 2x2 least-squares Jacobian per vertex, fitted on its 1-ring edges). This is the nsp atlas's
% canonical frame (nsp/atlas/frames.py, domain "sphere") in closed form: on the sphere the trivial
% connection with +1 at both poles IS the meridian frame (rheome.operators.gauge: 0.36 deg median; frames.py
% checks <= 0.2 deg), so solving for it would only add solver error.
%
% ⭐ Use this, not rheome.operators.gauge's default poles, whenever frame COMPONENTS are compared across
% subjects: rheome.operators.gauge places its poles per mesh (and an older convention put them where alpha
% power was lowest), so their location -- and every component -- differs between subjects.
%
% ⚠ THE FRAME IS MEANINGLESS NEAR A POLE, AND +z IS IN ACTIVE CORTEX (precentral: mu/alpha). The
% meridians converge, so a patch of half-width rho at colatitude theta sees the frame turn by about
% 2*rho*cot(theta) relative to parallel transport; a pooled vector there mixes directions. Two
% remedies, both measured by rheome.geom.spherepatches: drop the patch (.excluded) or read it in the
% second chart, Axis=[1 0 0], whose poles lie on the equator of the first. The two charts differ
% by a known per-vertex rotation (.e1 of both), so nothing is lost by switching.
%
% INPUTS:
%   V [nV x 3] cortex vertices, F [nF x 3] faces, N [nV x 3] cortex vertex normals (outward)
%   Sph [nV x 3] the registered-sphere point of each vertex (any radius; S.Sphere(gv,:))
%   Axis           pole axis on the sphere, default [0 0 1]
%   MaxCondition   a vertex whose push-forward Jacobian has condition above this is singular (50,
%                  as frames.py)
%
% OUTPUT (struct fr):
%   .e1 .e2 [nV x 3]  north and west (e2 = N x e1), unit, tangent; NaN at singular vertices
%   .colat  [nV x 1]  colatitude on the sphere from +Axis, radians    .lon [nV x 1] longitude, rad
%   .north .south     the pole vertices (closest to +Axis / -Axis)
%   .singular [nV x 1] logical: the poles and ill-conditioned push-forwards
%   .condition [nV x 1]  .axis
%
% See also: rheome.geom.spherepatches, rheome.operators.gauge, rheome.operators.trivial_connection, rheome.scale.measure_atlas
%
% Author: Diellor Basha, 2026

    arguments
        V (:,3) double
        F (:,3) double
        N (:,3) double
        Sph (:,3) double
        opts.Axis (1,3) double = [0 0 1]
        opts.MaxCondition (1,1) double = 50
    end
    nV = size(V,1);
    ax = opts.Axis / norm(opts.Axis);
    u = Sph ./ vecnorm(Sph,2,2);
    c = u*ax';
    mer = ax - c.*u;                                     % the meridian towards +Axis, tangent
    es = mer ./ max(vecnorm(mer,2,2), eps);
    [a1, a2] = i_basis(ax);
    lon = atan2(u*a2', u*a1');

    % per-vertex Jacobian sphere -> cortex, least squares over the directed 1-ring edges
    [sb1, sb2] = i_tangent(u);  [db1, db2] = i_tangent(N ./ vecnorm(N,2,2));
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    A = [E(:,1); E(:,2)];  B = [E(:,2); E(:,1)];
    ds = u(B,:) - u(A,:);  dd = V(B,:) - V(A,:);
    s1 = sum(ds.*sb1(A,:),2);  s2 = sum(ds.*sb2(A,:),2);
    c1 = sum(dd.*db1(A,:),2);  c2 = sum(dd.*db2(A,:),2);
    acc = @(x) accumarray(A, x, [nV 1]);
    S11 = acc(s1.*s1); S12 = acc(s1.*s2); S22 = acc(s2.*s2);
    C11 = acc(c1.*s1); C12 = acc(c1.*s2); C21 = acc(c2.*s1); C22 = acc(c2.*s2);
    dS = S11.*S22 - S12.^2;                              % inv(SS) = [S22 -S12; -S12 S11]/dS
    J11 = (C11.*S22 - C12.*S12)./dS;  J12 = (C12.*S11 - C11.*S12)./dS;
    J21 = (C21.*S22 - C22.*S12)./dS;  J22 = (C22.*S11 - C21.*S12)./dS;
    m1 = sum(es.*sb1,2);  m2 = sum(es.*sb2,2);
    e1 = (J11.*m1 + J12.*m2).*db1 + (J21.*m1 + J22.*m2).*db2;
    e1 = e1 ./ max(vecnorm(e1,2,2), eps);
    % condition of each 2x2 J from its singular values
    tr = J11.^2 + J12.^2 + J21.^2 + J22.^2;  dt = abs(J11.*J22 - J12.*J21);
    disc = sqrt(max(tr.^2 - 4*dt.^2, 0));
    cond = sqrt((tr + disc) ./ max(tr - disc, realmin));

    [~, north] = max(c);  [~, south] = min(c);
    sing = ~(cond <= opts.MaxCondition);  sing([north south]) = true;
    e1(sing,:) = NaN;
    fr = struct('e1', e1, 'e2', cross(N ./ vecnorm(N,2,2), e1, 2), 'colat', acos(max(-1,min(1,c))), ...
                'lon', lon, 'north', north, 'south', south, 'singular', sing, 'condition', cond, 'axis', ax);
end

function [b1, b2] = i_tangent(n)
% any orthonormal pair spanning each tangent plane (frames.py tangent_basis)
    ref = repmat([1 0 0], size(n,1), 1);  ref(abs(n(:,1)) >= 0.9, :) = repmat([0 1 0], sum(abs(n(:,1)) >= 0.9), 1);
    b1 = ref - sum(ref.*n,2).*n;  b1 = b1 ./ vecnorm(b1,2,2);  b2 = cross(n, b1, 2);
end

function [a1, a2] = i_basis(ax)
% the longitude reference: for +z it is (x, y), so lon = atan2(y, x) as on fsaverage
    if abs(ax(3)) > 0.9, a1 = [1 0 0]; else, a1 = [0 0 1]; end
    a1 = a1 - (a1*ax')*ax;  a1 = a1/norm(a1);  a2 = cross(ax, a1);
end

% Author: Diellor Basha, 2026
