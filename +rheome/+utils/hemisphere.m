function [Sh, dbh] = hemisphere(S, dbG, h)
% UTILS.HEMISPHERE  Restrict a cortex surface + Dirac basis to ONE hemisphere.
%
%   [Sh, dbh] = rheome.utils.hemisphere(S, dbG, h)
%   Sh        = rheome.utils.hemisphere(S, [], h)      % surface only
%
% The relative Dirac is built per hemisphere, so per-hemisphere work is the
% natural unit. This returns the chosen hemisphere as its OWN surface + basis, re-indexed
% 1..nVh, so downstream seeding / filtering / show.field need no masking.
%
% INPUTS:
%   S   : full-cortex surface (rheome.io.read.surface) -- must carry S.Hemi (Structures-atlas split)
%   dbG : full relative-Dirac eigenbasis (rheome.eigen.dirac_frame), or [] for surface only
%   h   : hemisphere index (1 or 2) OR label ('L' / 'R')
%
% OUTPUTS:
%   Sh  : hemisphere surface (.Vertices/.Faces re-indexed, .VertNormals sliced, .nV/.nF,
%         .GlobalVertices = local->global map, .Comment tagged with the hemisphere)
%   dbh : Dirac basis restricted to that hemisphere (local vertex indices), B-orthonormal
%
% See also: rheome.io.read.surface, rheome.eigen.dirac_frame, show.field
%
% Author: Diellor Basha, 2026

    if ischar(h) || (isstring(h) && isscalar(h))
        h = find(strcmpi(S.HemiLabel, sprintf('Cortex %s', upper(char(h)))), 1);
        if isempty(h), error('utils:hemisphere:label', 'No hemisphere with that label in S.HemiLabel.'); end
    end
    if ~isfield(S,'Hemi') || isempty(S.Hemi) || h < 1 || h > numel(S.Hemi)
        error('utils:hemisphere:noSplit', ...
            'Surface has no hemisphere split. Read it via rheome.io.read.surface (needs the Structures atlas).');
    end
    vidx  = S.Hemi{h}(:);
    remap = zeros(S.nV, 1);  remap(vidx) = 1:numel(vidx);

    % --- sub-surface (faces fully inside this hemisphere, re-indexed 1..nVh) ---
    keep = all(ismember(S.Faces, vidx), 2);
    Sh = struct();
    Sh.Vertices    = S.Vertices(vidx, :);
    Sh.Faces       = remap(S.Faces(keep, :));
    Sh.nV          = numel(vidx);
    Sh.nF          = size(Sh.Faces, 1);
    Sh.VertNormals = [];  if ~isempty(S.VertNormals), Sh.VertNormals = S.VertNormals(vidx, :); end
    Sh.Comment     = sprintf('%s | %s', S.Comment, S.HemiLabel{h});
    Sh.Hemi = {};  Sh.HemiLabel = {};                     % already a single hemisphere
    Sh.SurfaceFile    = S.SurfaceFile;
    Sh.GlobalVertices = vidx;                             % local -> global vertex map

    % --- sub-basis: this hemisphere's modes on this hemisphere's rows ---
    dbh = [];
    if nargin >= 2 && ~isempty(dbG)
        cols = find(dbG.Hemisphere == h);                 % modes belonging to hemisphere h
        rows = reshape((vidx' - 1)*4 + (1:4)', [], 1);    % its 4-quaternion vertex blocks (in vidx order)
        dbh = struct();
        dbh.Phi        = dbG.Phi(rows, cols);             % [4nVh x mh], local vertex order
        dbh.Mass       = dbG.Mass(rows, rows);
        dbh.Lambda     = dbG.Lambda(cols);
        dbh.nVert      = numel(vidx);
        dbh.nModes     = numel(cols);
        dbh.Tau        = dbG.Tau;
        dbh.Hemisphere = ones(numel(cols), 1);
        if isfield(dbG,'Scales'),    dbh.Scales    = dbG.Scales(h, :); end   % this hemisphere's [sL sE]
        if isfield(dbG,'Normalize'), dbh.Normalize = dbG.Normalize;    end
    end
end

% Author: Diellor Basha, 2026
