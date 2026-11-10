% REGISTER  A subject registered to its template (fsaverage): the one operation of subject to template, in two bases.
%
% The registration is the sphere registration R [template vertices x subject vertices] (the two registered spheres). It is
% applied to something measured on the subject's own cortex (never re-derived on the template's smaller one), by VERTICES (a
% map, points) or in the MODES (R written in the two Hodge-Dirac bases, T = A_template R P_subject: a map, or the weights of a
% kernel). The same package, names and mathematics as nxr.register (nxr) and the Cortical Flow Editor's "Register to template".
%
% The sphere registration (nxr.register.sphere):
%   rheome.register.hemispheres          - 'L'/'R' per vertex: the two largest pieces, told apart by the left axis (+y SCS, -x MNI)
%   rheome.register.registration         - R: points on a registered sphere sampled from a source's cortex (barycentric, nearest, shepard)
%   rheome.register.shepard_weights      - Brainstorm's Shepard weights of points from points
%
% By vertices (nxr.register.vertices):
%   rheome.register.register_values      - a scalar map by R, NaN where a template point sees no subject vertex
%   rheome.register.nearest_vertices     - points to their nearest template vertex on the registered sphere
%
% In the modes (nxr.register.modes):
%   rheome.register.solve                - one Hodge-Dirac solve with its glue (nxr's Solve)
%   rheome.register.surfaces             - the two surfaces, spheres and sampling a registration is built from (nxr's Surfaces)
%   rheome.register.glued_registration   - R on the two glued surfaces (bridge rows from the cut walls)
%   rheome.register.mode_registration    - T [K_template x K_subject] (and R): R written in the two mode bases
%   rheome.register.check_mode_registration - T measured on smooth fields (spherical harmonics l = 1..3)
%   rheome.register.register_in_modes    - a scalar map registered through T
%   rheome.register.to_glued, to_input   - values between the input vertices and the glued surface
%
% The stored registration of a dataset is read by rheome.io.read.nxrregistration; a subject's fused weights are registered by
% rheome.inverse.hodge_register. Not ported: nxr.register.tangent (tangent fields) and nxr.register.interpolation.
%
% Author: Diellor Basha, 2026
