function c = cortexmask(src, nV)
% GEOM.CORTEXMASK  Which vertices are cortex: the medial wall closes the hemisphere but holds no grey matter.
%
%   c = rheome.geom.cortexmask(name, nV) % from the cached Desikan-Killiany atlas (+data/<name>/atlas.mat)
%   c = rheome.geom.cortexmask(A, nV)    % from one atlas struct (.Label, .Membership [nParcel x nV])
%   c = rheome.geom.cortexmask(S)        % the mask a surface carries (S.Cortex), else all true
%
% FreeSurfer closes each hemisphere across the corpus callosum and the diencephalic cut so it can be
% inflated to a sphere; that patch is not cortex (thickness 0, no Desikan-Killiany label). The source
% space still puts dipoles on it, so values there are leakage, not measurements. Rheome keeps the
% closed mesh for its operators (rheome.operators.harmonicfill) and EXCLUDES the wall from every
% measurement: areas, sums, means, counts, plants and spin nulls run over c only (decision 3882c4ce,
% Diellor 2026-10-14; review: medial-wall-review.md).
%
% Cortex = the vertices some parcel covers, other than 'unknown', 'corpuscallosum' and 'Medial_wall'.
% The name form reads Desikan-Killiany, falling back to Destrieux; with neither cached it returns all
% true(nV, 1) and warns (rheome:cortexmask:none), so a surface without an atlas behaves as before.
%
% ⚠ Labels are matched by prefix, not substring: 'medialorbitofrontal' is cortex. The older filters in
%   measure_atlas / template matched "medial" anywhere and dropped it (fixed with this function).
%
% OUTPUT
%   c   [nV x 1] logical, true on cortex
%
% See also: rheome.operators.harmonicfill, rheome.geom.tiles, rheome.load.surface, rheome.load.bases
%
% Author: Diellor Basha, 2026

    if isstruct(src) && isfield(src, 'Vertices')                % a surface
        nv = size(src.Vertices, 1);
        if isfield(src, 'Cortex') && numel(src.Cortex) == nv, c = logical(src.Cortex(:)); else, c = true(nv, 1); end
        return
    end
    if ischar(src) || isstring(src)                            % a cached dataset
        name = char(src);  A = [];
        for at = ["Desikan-Killiany" "Destrieux"]
            try, A = rheome.load.atlas(name, at); break; catch, end %#ok<CTCH>
        end
        if isempty(A)
            warning('rheome:cortexmask:none', '%s: no Desikan-Killiany or Destrieux atlas cached; the medial wall is not masked.', name);
            c = true(nV, 1);  return
        end
        src = A;
    end
    lab = lower(strtrim(string(src.Label(:))));
    wall = startsWith(lab, ["unknown" "corpuscallosum" "corpus callosum" "medial_wall" "medial wall"]);
    c = full(any(src.Membership(~wall, :), 1)).';
    if numel(c) < nV, c(end+1:nV) = false; end                 % parcels stored on fewer columns
    c = c(1:nV);
end

% Author: Diellor Basha, 2026
