function A = atlas(SurfaceFile)
% IO.READ.ATLAS  Read every parcellation (atlas) from a Brainstorm cortical surface.
%
%   A = rheome.io.read.atlas(SurfaceFile)
%
% Returns a [1 x nAtlas] struct array, one entry per atlas stored in the surface. This is
% the ROI AXIS of the cortical-flow tensor: `Membership` is the reduction operator, so
% averaging a vertex field over ROIs is one sparse multiply rather than a loop.
%
% OUTPUT (per atlas):
%   .Name        atlas name ('Desikan-Killiany', 'Destrieux', 'Structures', ...)
%   .nScout      number of scouts
%   .Label       {1 x nScout} scout labels
%   .Region      {1 x nScout} Brainstorm region codes ('LT','RF','LO', ...)
%   .Hemi        {1 x nScout} 'L' | 'R', the first character of .Region
%   .Vertices    {1 x nScout} COLUMN vectors of 1-based vertex indices
%   .Seed        [1 x nScout] seed vertex of each scout
%   .nVertex     [1 x nScout] vertices per scout
%   .Membership  sparse logical [nScout x nV] -- scout s contains vertex v
%   .nCovered    vertices belonging to at least one scout
%   .Coverage    nCovered / nV
%   .Disjoint    true if no vertex belongs to two scouts
%   .nV          vertices in the surface
%
% ⚠ SCOUTS CARRY `Label`, NOT `Name`. Unlike almost every other Brainstorm struct. Reading
% `.Name` errors with "Unrecognized field name" -- a fine failure, but only if you were
% looking. (`.Name` is accepted as a fallback if a file ever has it.)
%
% ⚠ VERTICES ARE STORED AS ROW VECTORS. `vertcat(sc.Vertices)` over a scout array errors
% ("Dimensions of arrays being concatenated are not consistent") the moment orientation is
% mixed. Everything here is normalised to COLUMNS.
%
% ⚠ A PARCELLATION IS NOT A PARTITION. Desikan-Killiany covers 18743 of 20484 vertices on
% these surfaces; the medial wall belongs to no scout. Summing over ROIs is therefore NOT
% the whole-cortex total, and nothing announces the difference -- check `.Coverage`.
% Overlap is likewise possible (`.Disjoint`), in which case ROI sums double-count.
%
% ⚠ AREA WEIGHTING LIVES OUTSIDE THIS FUNCTION. An atlas is anatomy; the weights are the
% surface's. Compose them at the point of use:
%     w     = full(diag(M));                     % lumped vertex areas
%     P     = double(A.Membership);              % [nScout x nV]
%     areas = P * w;
%     roiMean = (P * (w .* x)) ./ areas;         % area-weighted ROI mean of x
%
% See also: rheome.import.atlas, rheome.load.atlas, rheome.io.read.surface
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(SurfaceFile)
        error('io:read:atlas:args', 'A surface file path is required.');
    end
    if ~exist(SurfaceFile, 'file')
        error('io:read:atlas:notFound', 'Surface file not found: %s', SurfaceFile);
    end

    raw = builtin('load', SurfaceFile, 'Atlas', 'Vertices');
    if ~isfield(raw, 'Atlas') || isempty(raw.Atlas)
        error('io:read:atlas:noAtlas', ...
            'No Atlas in %s -- the surface carries no parcellation.', SurfaceFile);
    end
    if ~isfield(raw, 'Vertices') || isempty(raw.Vertices)
        error('io:read:atlas:noVertices', ...
            'No Vertices in %s; cannot size the membership matrix.', SurfaceFile);
    end
    nV = size(raw.Vertices, 1);

    R = raw.Atlas;
    A = repmat(i_empty(nV), 1, numel(R));
    for k = 1:numel(R)
        A(k) = i_one(R(k), nV, SurfaceFile);
    end
end

% ---- one atlas ----
function a = i_one(rec, nV, srcFile)
    a      = i_empty(nV);
    a.Name = i_get(rec, 'Name', '');
    sc     = i_get(rec, 'Scouts', []);
    nS     = numel(sc);
    a.nScout = nS;
    if nS == 0, return; end

    a.Label    = cell(1, nS);   a.Region = cell(1, nS);
    a.Hemi     = cell(1, nS);   a.Vertices = cell(1, nS);
    a.Seed     = zeros(1, nS);  a.nVertex  = zeros(1, nS);

    rows = cell(1, nS);  cols = cell(1, nS);
    for j = 1:nS
        % Label, not Name -- see the header. Name is only a fallback.
        lab = i_get(sc(j), 'Label', '');
        if isempty(lab), lab = i_get(sc(j), 'Name', sprintf('scout%d', j)); end
        a.Label{j} = char(lab);

        reg = char(i_get(sc(j), 'Region', ''));
        a.Region{j} = reg;
        if isempty(reg), a.Hemi{j} = ''; else, a.Hemi{j} = upper(reg(1)); end

        v = double(i_get(sc(j), 'Vertices', []));
        v = v(:);                                        % ROW -> COLUMN, always
        if any(v < 1 | v > nV | mod(v, 1) ~= 0)
            error('io:read:atlas:vertexIndex', ...
                'Atlas ''%s'' scout ''%s'' has out-of-range vertex indices (nV = %d) in %s.', ...
                a.Name, a.Label{j}, nV, srcFile);
        end
        a.Vertices{j} = v;
        a.nVertex(j)  = numel(v);
        a.Seed(j)     = i_get(sc(j), 'Seed', NaN);

        rows{j} = repmat(j, numel(v), 1);
        cols{j} = v;
    end

    r = vertcat(rows{:});  c = vertcat(cols{:});
    a.Membership = sparse(r, c, true, nS, nV);
    a.nCovered   = numel(unique(c));
    a.Coverage   = a.nCovered / nV;
    % sparse() SUMS duplicates, so a vertex in two scouts appears twice in c but once as a
    % column of Membership -- compare the raw count, not the matrix's nnz.
    a.Disjoint   = (numel(c) == a.nCovered);
end

function a = i_empty(nV)
    a = struct('Name', '', 'nScout', 0, 'Label', {{}}, 'Region', {{}}, 'Hemi', {{}}, ...
               'Vertices', {{}}, 'Seed', [], 'nVertex', [], ...
               'Membership', sparse(false(0, nV)), 'nCovered', 0, 'Coverage', 0, ...
               'Disjoint', true, 'nV', nV);
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
