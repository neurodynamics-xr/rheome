function [order, t, wind] = spiralorder(Sph, opts)
% GEOM.SPIRALORDER  A 1-D order of sphere points along a pole-to-pole spiral: vertices or tiles.
%
%   [order, t] = rheome.geom.spiralorder(Sph)                         % 16 turns, +z to -z, even in colatitude
%   [order, t] = rheome.geom.spiralorder(Sph, Turns=32, Kind="area")
%   order = rheome.geom.spiralorder(P.centre)                         % sort rheome.geom.spherepatches tiles
%
% The spiral starts at the gauge's north pole (+Axis, rheome.geom.sphereframe's .north) and winds Turns
% times eastward (longitude increasing, as rheome.geom.sphereframe's .lon) to the south pole. Each point
% is projected onto it along its own meridian, onto the nearest winding, and t is where it lands:
% Sph(order,:) runs down the spiral. Because t is a function of the REGISTERED sphere point only,
% a given t is the same anatomical place in every subject (Sph = S.Sphere(gv,:), sphere.reg).
%
% Kind sets how the windings are spaced, g(colatitude) in [0, 1] growing linearly along the spiral:
%   "colat"      g = theta/pi: a constant gap of 180/Turns degrees between windings, so a winding
%                is a strip of constant width -- the locality choice (the default)
%   "area"       g = (1 - cos theta)/2: equal area per unit t, the spherical Fibonacci spiral's
%                spacing; windings thin out (in degrees) towards the equator
%   "loxodrome"  g linear in the Mercator ordinate, so the spiral keeps a constant bearing; it
%                winds infinitely often at a pole, so g is clamped at Cap degrees from each pole
%
% ⚠ A WINDING IS A STRIP, SO AN ORDER CAN'T BE LOCAL BOTH WAYS. Points one winding apart are
% neighbours on the sphere and about n/Turns apart in the order; that, not the curve, is the cost
% of any 1-D order of a 2-D sheet. Choose Turns for the use: few turns keep east-west neighbours
% close in the order, many turns keep each strip narrow.
%
% INPUTS:  Sph [n x 3] points on (or radially projectable to) the sphere, any radius
%   Turns  windings from pole to pole (positive, need not be an integer), default 16
%   Axis   pole axis, default [0 0 1] (the group gauge; [1 0 0] is its second chart)
%   Kind   "colat" (default) | "area" | "loxodrome"
%   Cap    degrees, loxodrome only: the polar cap that is flattened onto t = 0 or 1, default 2
%
% OUTPUT:
%   order [n x 1]  permutation: Sph(order,:) runs from the north pole to the south
%   t     [n x 1]  position on the spiral, 0 at north and 1 at south (a point within half a winding
%                  of a pole may land just outside [0, 1]; only its order matters)
%   wind  [n x 1]  the winding each point was projected onto, -1 .. ceil(Turns) (-1 only next to the north pole)
%
% See also: rheome.geom.sphereframe, rheome.geom.spherepatches, rheome.geom.tree
%
% Author: Diellor Basha, 2026

    arguments
        Sph (:,3) double
        opts.Turns (1,1) double {mustBePositive} = 16
        opts.Axis (1,3) double = [0 0 1]
        opts.Kind (1,1) string {mustBeMember(opts.Kind, ["colat" "area" "loxodrome"])} = "colat"
        opts.Cap (1,1) double {mustBePositive} = 2
    end
    ax = opts.Axis / norm(opts.Axis);
    u = Sph ./ vecnorm(Sph,2,2);
    c = max(-1, min(1, u*ax'));
    if abs(ax(3)) > 0.9, a1 = [1 0 0]; else, a1 = [0 0 1]; end   % as rheome.geom.sphereframe's longitude
    a1 = a1 - (a1*ax')*ax;  a1 = a1/norm(a1);  a2 = cross(ax, a1);
    f = mod(atan2(u*a2', u*a1'), 2*pi) / (2*pi);                 % longitude as a fraction of a turn
    switch opts.Kind
        case "colat", g = acos(c) / pi;
        case "area",  g = (1 - c) / 2;
        case "loxodrome"
            y0 = asinh(cotd(opts.Cap));
            g = (y0 - max(-y0, min(y0, asinh(c ./ sqrt(max(1 - c.^2, eps)))))) / (2*y0);
    end
    wind = round(opts.Turns*g - f);       % the winding passing closest along the meridian
    t = (wind + f) / opts.Turns;
    [~, order] = sort(t);
end

% Author: Diellor Basha, 2026
