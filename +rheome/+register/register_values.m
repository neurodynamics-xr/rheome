function y = register_values(R, x)
% REGISTER.REGISTER_VALUES  A subject's scalar maps registered by vertices: R applied, NaN where a template point has no subject vertex.
%
%   y = rheome.register.register_values(R, x)       % x [V_subject x m] measured on the subject's own cortex; y [V_template x m]
%
% A measurement is made on the subject's own cortex (physical units, the subject's lengths) and then registered, never the
% current itself re-derived on the template's smaller cortex (a size bias). NaN in x counts as 0. As nxr.register.vertices.register_values.
%
% See also: rheome.register.registration, rheome.register.nearest_vertices
%
% Author: Diellor Basha, 2026
arguments, R (:,:) double, x (:,:) double, end
x(isnan(x)) = 0;
y = R * x;
y(sum(R ~= 0, 2) == 0, :) = NaN;
end
