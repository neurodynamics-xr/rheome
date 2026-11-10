function reg = surfaces(subjectPositions, subjectFaces, subjectSphere, templatePositions, templateFaces, templateSphere, opts)
% REGISTER.SURFACES  What a registration is built from beyond the two solves: the two surfaces and their registered spheres.
%
%   reg = rheome.register.surfaces(subjectPositions, subjectFaces, subjectSphere, templatePositions, templateFaces, templateSphere)
%   reg = rheome.register.surfaces(..., SubjectLeft="+y", TemplateLeft="-x", Sampling="shepard", Neighbours=8, Power=2)
%
% The positions, faces (1-based) and registered spheres are of the INPUT surfaces (before the glue). SubjectLeft / TemplateLeft:
% the axis the left hemisphere lies on (rheome.register.hemispheres). Sampling: "barycentric", "nearest" or "shepard"
% (rheome.register.registration). reg is a struct with the fields subject_positions, ..., template_sphere, subject_left,
% template_left, sampling, neighbours, power. As nxr.register.modes.Surfaces.
%
% See also: rheome.register.solve, rheome.register.mode_registration
%
% Author: Diellor Basha, 2026
arguments
    subjectPositions (:,3) double, subjectFaces (:,3) double, subjectSphere (:,3) double
    templatePositions (:,3) double, templateFaces (:,3) double, templateSphere (:,3) double
    opts.SubjectLeft (1,1) string = "+y"
    opts.TemplateLeft (1,1) string = "-x"
    opts.Sampling (1,1) string {mustBeMember(opts.Sampling, ["barycentric" "nearest" "shepard"])} = "barycentric"
    opts.Neighbours (1,1) double = 8
    opts.Power (1,1) double = 2
end
reg = struct('subject_positions', subjectPositions, 'subject_faces', subjectFaces, 'subject_sphere', subjectSphere, ...
    'template_positions', templatePositions, 'template_faces', templateFaces, 'template_sphere', templateSphere, ...
    'subject_left', opts.SubjectLeft, 'template_left', opts.TemplateLeft, 'sampling', opts.Sampling, ...
    'neighbours', opts.Neighbours, 'power', opts.Power);
end
