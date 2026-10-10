function build_hcp_template(cortexFile, fibersFile, outFile)
% CONNECTOME.BUILD_HCP_TEMPLATE  One-time extraction of the HCP-1065 template into hcp_template.mat.
%
%   rheome.connectome.build_hcp_template(cortexFile, fibersFile, outFile)
%
% The HCP-1065 template is NOT part of stock Brainstorm -- it was manually added to a
% @default_subject. This script reads that (Brainstorm) template cortex + fibers ONCE and writes
% a self-contained hcp_template.mat, so rheome.connectome.resolve needs no Brainstorm at runtime.
% Reuses rheome.io.read.surface / rheome.io.read.fibers for a consistent read. The file is derived
% from WU-Minn HCP data: share it only under the HCP Open Access Data Use Terms.
%
% INPUTS:
%   cortexFile  a @default_subject cortex .mat WITH a Reg.Sphere, whose surface the fiber endpoints
%               sit on (tess_cortex_pial_low is closest for HCP-1065).
%   fibersFile  the sibling tess_fibers_*hcp1065*.mat
%   outFile     output hcp_template.mat (put it at fullfile(rheome.load.root(), 'hcp_template.mat'))
%
% Author: Diellor Basha, 2026

    templateSurface        = rheome.io.read.surface(cortexFile);
    templateVertices       = templateSurface.Vertices;              % [numTemplateVerts x 3]
    templateSphereVertices = templateSurface.Sphere;               % [numTemplateVerts x 3] FreeSurfer sphere
    templateHemi           = templateSurface.Hemi;                 % {leftVerts, rightVerts}
    assert(~isempty(templateSphereVertices), 'Template cortex has no Reg.Sphere.');
    assert(numel(templateHemi) == 2, 'Template cortex has no left/right hemisphere split.');

    templateFiberEndpoints = rheome.io.read.fibers(fibersFile);           % [numFibers x 2 x 3]

    [~, cn, ce] = fileparts(cortexFile);  [~, fn, fe] = fileparts(fibersFile);   % names only: no local paths
    provenance = struct('cortexFile', [cn ce], 'fibersFile', [fn fe], ...
        'numFibers', size(templateFiberEndpoints,1), 'numTemplateVerts', size(templateVertices,1));

    save(outFile, 'templateVertices', 'templateSphereVertices', 'templateHemi', ...
        'templateFiberEndpoints', 'provenance', '-v7');
    fprintf('wrote %s: %d verts, %d fibers, hemis [%d %d]\n', outFile, ...
        size(templateVertices,1), size(templateFiberEndpoints,1), numel(templateHemi{1}), numel(templateHemi{2}));
end

% Author: Diellor Basha, 2026
