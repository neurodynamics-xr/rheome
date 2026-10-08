function S = surface(SurfaceFile)
% IO.READ.SURFACE  Read a Brainstorm cortical surface .mat into a plain struct.
%
%   S = rheome.io.read.surface(SurfaceFile)
%
% A Brainstorm surface file is just a MATLAB .mat holding, among others, the
% variables 'Vertices' [nV x 3] (in METERS) and 'Faces' [nF x 3] (1-based vertex
% indices). This reader loads those, validates them, and returns a minimal,
% self-describing struct so the rest of the demo never has to know anything about
% the Brainstorm environment (no bst_get, no protocol database, no GUI).
%
% INPUT:
%   SurfaceFile : path to a Brainstorm surface .mat (e.g. tess_cortex_pial_low.mat)
%
% OUTPUT (struct S):
%   .Vertices    [nV x 3] vertex coordinates (meters)
%   .Faces       [nF x 3] triangle vertex indices (1-based)
%   .nV, .nF     counts
%   .VertNormals [nV x 3] vertex normals if present in the file, else []
%   .Sphere      [nV x 3] FreeSurfer registration-sphere coords (Reg.Sphere), else []
%   .VertConn    [nV x nV] sparse vertex adjacency if present, else []
%   .Comment     the surface's Brainstorm comment if present, else ''
%   .Hemi        {leftVerts, rightVerts} vertex-index lists from the 'Structures' atlas
%                (Cortex L / Cortex R scouts) -- the AUTHORITATIVE hemisphere split; {} if
%                the atlas is absent or does not partition the surface
%   .HemiLabel   {1 x 2} labels for .Hemi ({'Cortex L','Cortex R'})
%   .SurfaceFile the input path (provenance)
%
% See also: rheome.operators.mass, rheome.operators.laplace_beltrami, rheome.eigen.dirac_frame
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(SurfaceFile)
        error('io:read:surface:args', 'A surface file path is required.');
    end
    if ~exist(SurfaceFile, 'file')
        error('io:read:surface:notFound', 'Surface file not found: %s', SurfaceFile);
    end

    % A Brainstorm surface .mat stores its fields as top-level variables.
    raw = load(SurfaceFile);
    if ~isfield(raw, 'Vertices') || ~isfield(raw, 'Faces')
        error('io:read:surface:badFile', ...
            'File does not look like a Brainstorm surface (missing Vertices/Faces): %s', SurfaceFile);
    end

    V = double(raw.Vertices);
    F = double(raw.Faces);

    % --- validate geometry ---
    if size(V, 2) ~= 3
        error('io:read:surface:vertices', 'Vertices must be [nV x 3], got [%d x %d].', size(V,1), size(V,2));
    end
    if size(F, 2) ~= 3
        error('io:read:surface:faces', 'Faces must be triangles [nF x 3], got [%d x %d].', size(F,1), size(F,2));
    end
    nV = size(V, 1);
    if any(F(:) < 1) || any(F(:) > nV) || any(mod(F(:), 1) ~= 0)
        error('io:read:surface:faceIndex', 'Faces contain out-of-range or non-integer vertex indices.');
    end

    S = struct();
    S.Vertices    = V;
    S.Faces       = F;
    S.nV          = nV;
    S.nF          = size(F, 1);
    S.VertNormals = [];
    if isfield(raw, 'VertNormals') && ~isempty(raw.VertNormals) && size(raw.VertNormals,1) == nV
        S.VertNormals = double(raw.VertNormals);
    end
    S.Sphere = [];
    if isfield(raw, 'Reg') && isstruct(raw.Reg) && isfield(raw.Reg, 'Sphere') ...
            && isfield(raw.Reg.Sphere, 'Vertices') && size(raw.Reg.Sphere.Vertices,1) == nV
        S.Sphere = double(raw.Reg.Sphere.Vertices);
    end
    S.VertConn = [];
    if isfield(raw, 'VertConn') && ~isempty(raw.VertConn) && size(raw.VertConn,1) == nV
        S.VertConn = raw.VertConn;
    end
    S.Comment = '';
    if isfield(raw, 'Comment'), S.Comment = raw.Comment; end

    % --- hemisphere vertex lists from the 'Structures' atlas (authoritative Cortex L/R
    %     scouts). This is the CORRECT split -- conncomp on the mesh is fragile
    %     (non-manifold edges / stray components can mis-split). Only set when the L/R
    %     scouts cleanly partition every vertex; otherwise leave {} (callers fall back).
    S.Hemi = {};  S.HemiLabel = {};
    if isfield(raw, 'Atlas') && ~isempty(raw.Atlas)
        iStr = find(strcmpi({raw.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStr)
            sc = raw.Atlas(iStr).Scouts;
            reg = repmat(' ', 1, numel(sc));                          % hemisphere char per scout
            for k = 1:numel(sc)
                if isfield(sc(k),'Region') && ~isempty(sc(k).Region), reg(k) = upper(sc(k).Region(1)); end
            end
            hemi = {}; lab = {};
            for h = 'LR'
                vv = [];
                for k = find(reg == h), vv = [vv; sc(k).Vertices(:)]; end %#ok<AGROW>
                if ~isempty(vv), hemi{end+1} = unique(vv); lab{end+1} = sprintf('Cortex %s', h); end %#ok<AGROW>
            end
            if numel(hemi) == 2 && numel(unique(vertcat(hemi{:}))) == nV   % clean full partition
                S.Hemi = hemi;  S.HemiLabel = lab;
            end
        end
    end

    S.SurfaceFile = SurfaceFile;
end

% Author: Diellor Basha, 2026
