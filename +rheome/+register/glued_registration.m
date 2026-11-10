function R = glued_registration(subject, template, reg)
% REGISTER.GLUED_REGISTRATION  The registration on the glued surfaces: R [V_template x V_subject] sparse.
%
%   R = rheome.register.glued_registration(subject, template, reg)       % solves: rheome.register.solve; reg: rheome.register.surfaces
%
% The template's kept vertices are pulled from the subject's cortex through the registered spheres (rheome.register.registration,
% a hemisphere at a time; the subject's bridge is not cortex), each template bridge vertex is the blend of its 3 nearest
% cut-out wall vertices, those pulled the same way. As nxr.register.modes.glued_registration.
%
% See also: rheome.register.mode_registration
%
% Author: Diellor Basha, 2026
arguments, subject (1,1) struct, template (1,1) struct, reg (1,1) struct, end
ms = subject.kept; mt = template.kept;
if max(abs(reg.template_positions(template.keep(1:mt), :) - template.positions(1:mt, :)), [], 'all') > 1e-6
    error('rheome:register:surface', 'the template solve is not on the template surface given (positions differ)');
end
if max(abs(reg.subject_positions(subject.keep(1:ms), :) - subject.positions(1:ms, :)), [], 'all') > 1e-6
    error('rheome:register:surface', 'the subject solve is not on the subject surface given (positions differ)');
end
hs = rheome.register.hemispheres(reg.subject_positions, reg.subject_faces, reg.subject_left);
ht = rheome.register.hemispheres(reg.template_positions, reg.template_faces, reg.template_left);
nIn = size(reg.subject_positions, 1);
cortex = false(nIn, 1); cortex(subject.keep(1:ms)) = true;
cols = zeros(nIn, 1); cols(subject.keep(1:ms)) = 1:ms;
pull = @(pts) rheome.register.registration(reg.subject_sphere, reg.subject_faces, hs, reg.template_sphere(pts, :), ht(pts), ...
    SrcCortex=cortex, SrcColumns=cols, Sampling=reg.sampling, Neighbours=reg.neighbours, Power=reg.power);
[cut, B] = bridge_rows(template, reg.template_positions);
R = [pull(template.keep(1:mt)); B * pull(cut)];
R(size(template.positions, 1), size(subject.positions, 1)) = 0;     % the subject's bridge columns: never sampled
R = R(1:size(template.positions, 1), 1:size(subject.positions, 1));
end
