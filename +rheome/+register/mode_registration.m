function [T, R] = mode_registration(subject, template, reg)
% REGISTER.MODE_REGISTRATION  The registration in the mode bases: a subject's Hodge-Dirac coefficients registered onto a template's.
%
%   [T, R] = rheome.register.mode_registration(subject, template, reg)
%
% subject, template  rheome.register.solve (the basis and its glue);  reg  rheome.register.surfaces
% T  [K_template x K_subject]: coefficients c on the subject (its field f = P c) to the template's: T c. Block-diagonal by kind
%    (exact -> exact, coexact -> coexact): T = A_template R P_subject, with A the template's analysis (exact 2 P' M; coexact
%    G^-1 P' M, G = P' M P on the coexact modes, which are not orthonormal). A field projection, not mode matching: mode signs
%    and rotations inside near-degenerate eigenspaces cancel, and T is exact for what the two bands share.
% R  [V_template x V_subject] sparse, the registration on the glued surfaces (rheome.register.glued_registration).
%
% As nxr.register.modes.mode_registration (the Cortical Flow Editor's registration to the template). Its check:
% rheome.register.check_mode_registration.
%
% See also: rheome.register.check_mode_registration, rheome.inverse.hodge_register, rheome.register.registration
%
% Author: Diellor Basha, 2026
arguments, subject (1,1) struct, template (1,1) struct, reg (1,1) struct, end
R = rheome.register.glued_registration(subject, template, reg);
T = fold(template, R * double(subject.P)) .* (template.kind(:) == subject.kind(:)');
end
