function op = operator(S, bases)
% DETECT.OPERATOR  Precompute the (mesh-only) operators the detectors reuse across frames.
%
%   op = rheome.detect.operator(S)             % build the operators from the mesh
%   op = rheome.detect.operator(S, bases)      % reuse the cached connection from rheome.load.bases (no rebuild)
%
% Every detect function needs the same J-INDEPENDENT operators: the whole-surface face gradient
% (for divergence / curl) and, per hemisphere, the Levi-Civita CONNECTION LAPLACIAN (tangent
% frames + parallel-transport rotations, for the winding number). rheome.detect.criticalPoints /
% rheome.detect.vortex rebuild these on every call; over a time series that repeats the SAME
% full-resolution factorisation for each frame. Build it ONCE here and pass op into the
% detectors -- e.g.
%       op = rheome.detect.operator(S, D.bases);
%       for t = 1:nT,  v(t) = rheome.detect.vortex(J(:,t), S, 'size', op);  end
%
% If a cached 'bases' (rheome.load.bases / D.bases) is supplied, each hemisphere's connection operator
% is taken from it (matched by global vertex set) instead of rebuilt -- the high-resolution,
% pre-computed operator. Missing / unmatched hemispheres fall back to a fresh build (identical
% result -- the connection Laplacian is deterministic, verified bit-for-bit against the cache).
%
% INPUT:
%   S      full-cortex OR single-hemisphere surface (rheome.io.read.surface / rheome.utils.hemisphere): needs
%          .Vertices, .Faces, .nV, .Hemi for the split (absent -> whole surface). A single-
%          hemisphere S should carry .GlobalVertices so it can be matched to the cache.
%   bases  (optional) rheome.load.bases(name) / D.bases -- per-hemisphere cached operators & eigenbases.
%
% OUTPUT (struct op):
%   .fg        rheome.operators.face_gradient(S) -- whole-surface gradient (divergence / curl)
%   .nH        number of hemispheres detected
%   .gv  {nH}  vertex indices of each hemisphere INTO the field J (global for a whole cortex,
%              1..nVh for a single-hemisphere S)
%   .Sh  {nH}  each hemisphere's surface (.Vertices/.Faces; centroids for the winding)
%   .C   {nH}  rheome.operators.connection_laplacian of each hemisphere (.e1 .e2 .Rt ...), cached or fresh
%   .nV        S.nV (guards J against a mismatched surface)
%   .cached [1 x nH] logical -- was this hemisphere's connection taken from the cache?
%
% See also: rheome.detect.criticalPoints, rheome.detect.vortex, rheome.load.bases, rheome.operators.connection_laplacian
%
% Author: Diellor Basha, 2026

    if nargin < 2, bases = []; end
    op = struct();
    op.nV = S.nV;
    op.fg = rheome.operators.face_gradient(S.Vertices, S.Faces);       % whole-surface, for div/curl

    if isfield(S,'Hemi') && ~isempty(S.Hemi), nH = numel(S.Hemi); else, nH = 1; end
    op.nH = nH;
    op.gv = cell(1, nH);  op.Sh = cell(1, nH);  op.C = cell(1, nH);  op.cached = false(1, nH);
    for h = 1:nH
        if nH == 1 && (~isfield(S,'Hemi') || isempty(S.Hemi))
            Sh = S;  gv = (1:S.nV)';
            gvGlobal = [];  if isfield(S,'GlobalVertices'), gvGlobal = S.GlobalVertices(:); end
        else
            Sh = rheome.utils.hemisphere(S, [], h);  gv = Sh.GlobalVertices;
            gvGlobal = gv;
        end
        op.Sh{h} = Sh;
        op.gv{h} = gv;

        C = i_cachedconn(bases, gvGlobal);                     % try the cache first
        if isempty(C)
            C = rheome.operators.connection_laplacian(Sh.Vertices, Sh.Faces);   % full-mesh, fresh
        else
            op.cached(h) = true;
        end
        op.C{h} = C;
    end
end

% ----- find a cached connection operator whose hemisphere matches this global vertex set -----
function C = i_cachedconn(bases, gvGlobal)
    C = [];
    if isempty(bases) || isempty(gvGlobal) || ~isstruct(bases), return; end
    for lab = {'L','R'}
        if isfield(bases, lab{1}) && isfield(bases.(lab{1}), 'gv') && isfield(bases.(lab{1}), 'conn')
            gvb = bases.(lab{1}).gv(:);
            if numel(gvb) == numel(gvGlobal) && isequal(gvb, gvGlobal)
                C = bases.(lab{1}).conn.C;  return;
            end
        end
    end
end

% Author: Diellor Basha, 2026
