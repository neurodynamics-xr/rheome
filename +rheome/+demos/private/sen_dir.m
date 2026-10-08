function d = sen_dir(pl)
% ⚠ OFF THE LATTICE AXIS, DELIBERATELY. A wave running along a grid axis gives every sensor
% in a column the SAME phase, so 64 sensors collapse to 8 distinct values and any per-sensor
% scatter shows 8 points. 30 degrees off-axis gives every sensor its own phase, which is both
% the more honest picture and the more general case.
    d = (cosd(30)*pl.e1 + sind(30)*pl.e2).';
    d = d / norm(d);
end

% Author: Diellor Basha, 2026
