function [subjectFiberEndpoints, provenance] = resolve(subjectName, subjectSurface)
% CONNECTOME.RESOLVE  Resolve the structural-connectome fiber source for a subject, with provenance.
%
%   [subjectFiberEndpoints, provenance] = rheome.connectome.resolve(subjectName, subjectSurface)
%
% Priority:
%   (1) the subject's OWN tractography if present  (+data/<name>/fibers_subject.mat)
%   (2) DEFAULT: the HCP-1065 group average (hcp_template.mat), registered to the subject
%       cortex through the FreeSurfer sphere -- for each fiber endpoint: nearest template-cortex
%       vertex -> its sphere coordinate -> nearest SUBJECT sphere vertex (per hemisphere).
%
% Returns the endpoints in SUBJECT coordinates [numFibers x 2 x 3] plus a provenance struct
% (fiberSource, registrationMethod, numFibers, templateSource) so the connectome basis is never a
% black box. Pure MATLAB; the default path needs only the template file, not Brainstorm.
%
% ⚠ The template is NOT part of the toolbox: it is derived from WU-Minn HCP data, which may only be
% redistributed under the HCP Open Access Data Use Terms, not under the toolbox's MIT licence. Download
% hcp_template.mat from the release page (https://github.com/neurodynamics-xr/rheome/releases), or build
% it from your own copy with rheome.connectome.build_hcp_template, and put it at
% fullfile(rheome.load.root(), 'hcp_template.mat') or point RHEOME_HCP_TEMPLATE at it.
%
% See also: rheome.connectome.build_hcp_template, rheome.operators.connectome, rheome.io.read.fibers
%
% Author: Diellor Basha, 2026

    subjectFiberFile = fullfile(rheome.load.root(), char(subjectName), 'fibers_subject.mat');
    if exist(subjectFiberFile, 'file')
        subjectFiberEndpoints = rheome.io.read.fibers(subjectFiberFile);
        provenance = struct('fiberSource','subject-specific', 'registrationMethod','none', ...
            'numFibers', size(subjectFiberEndpoints,1), 'templateSource','');
        return;
    end

    % ----- DEFAULT: HCP-1065 average, registered via the FreeSurfer sphere -----
    templateFile = getenv('RHEOME_HCP_TEMPLATE');
    if isempty(templateFile), templateFile = fullfile(rheome.load.root(), 'hcp_template.mat'); end
    if ~isfile(templateFile)
        error('connectome:notemplate', ['No HCP-1065 template at %s. It is distributed separately, under the ' ...
            'HCP Open Access Data Use Terms: download hcp_template.mat from ' ...
            'https://github.com/neurodynamics-xr/rheome/releases (or build it with ' ...
            'rheome.connectome.build_hcp_template) and put it there, or set RHEOME_HCP_TEMPLATE.'], templateFile);
    end
    template = load(templateFile);
    assert(~isempty(subjectSurface.Sphere), 'Subject has no Reg.Sphere; cannot register the HCP template.');

    numFibers = size(template.templateFiberEndpoints, 1);
    % stack the two ends explicitly: rows 1..numFibers = end 1, rows numFibers+1..2*numFibers = end 2
    firstEndXYZ  = squeeze(template.templateFiberEndpoints(:, 1, :));    % [numFibers x 3]
    secondEndXYZ = squeeze(template.templateFiberEndpoints(:, 2, :));    % [numFibers x 3]
    endpointXYZ  = [firstEndXYZ; secondEndXYZ];                          % [2*numFibers x 3]

    % 1) each endpoint -> nearest TEMPLATE cortex vertex -> its sphere coordinate
    nearestTemplateVertex = knnsearch(template.templateVertices, endpointXYZ);
    endpointSphereXYZ     = template.templateSphereVertices(nearestTemplateVertex, :);
    % 2) hemisphere of each endpoint (from the template hemi lists)
    isLeftEndpoint = ismember(nearestTemplateVertex, template.templateHemi{1});
    % 3) nearest SUBJECT sphere vertex within the SAME hemisphere
    subjectVertexIndex = zeros(size(endpointSphereXYZ,1), 1);
    for hemiSide = 1:2
        subjectHemiVertices = subjectSurface.Hemi{hemiSide};
        pick = (isLeftEndpoint == (hemiSide == 1));
        withinHemi = knnsearch(subjectSurface.Sphere(subjectHemiVertices,:), endpointSphereXYZ(pick,:));
        subjectVertexIndex(pick) = subjectHemiVertices(withinHemi);
    end
    % 4) subject-cortex coordinates, unstack back to [numFibers x 2 x 3]
    subjectEndpointXYZ = subjectSurface.Vertices(subjectVertexIndex, :);
    subjectFiberEndpoints = zeros(numFibers, 2, 3);
    subjectFiberEndpoints(:, 1, :) = subjectEndpointXYZ(1:numFibers, :);
    subjectFiberEndpoints(:, 2, :) = subjectEndpointXYZ(numFibers+1:end, :);

    provenance = struct('fiberSource','hcp-average', 'registrationMethod','freesurfer-sphere-nn', ...
        'numFibers', numFibers, 'templateSource', template.provenance.fibersFile);
end

% Author: Diellor Basha, 2026
