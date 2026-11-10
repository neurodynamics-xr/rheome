function z = register_in_modes(T, subject, template, x)
% REGISTER.REGISTER_IN_MODES  A scalar map on the subject's glued surface registered through T, in the modes.
%
%   z = rheome.register.register_in_modes(T, subject, template, x)      % x [V_subject_glued x m]; z [V_template_glued x m]
%
% x is analysed in the subject's modes, T applied, the result synthesised on the template. As nxr.register.modes.register_in_modes.
%
% See also: rheome.register.mode_registration, rheome.register.register_values
%
% Author: Diellor Basha, 2026
arguments, T (:,:) double, subject (1,1) struct, template (1,1) struct, x (:,:) double, end
z = double(template.P) * (T * fold(subject, x));
end
