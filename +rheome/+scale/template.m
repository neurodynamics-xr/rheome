function Tp = template(dmax, opts)
% SCALE.TEMPLATE  The shared template cortex and its group tiles: the space every participant is read in.
%
%   Tp = rheome.scale.template(3)                                  % RHEOME_TEMPLATE or RHEOME_BRAINSTORM's ICBM152
%   Tp = rheome.scale.template(4, Template='/path/ICBM152', Hemis=["L" "R"])
%
% A Brainstorm template folder gives tess_cortex_pial_low.mat (Reg.Sphere, Desikan-Killiany) and
% subjectimage_T1.mat (SCS, NCS). The group tiles are rheome.geom.tree on each template hemisphere, full to
% depth dmax; rheome.scale.grouptiles carries them to a participant through the registered sphere.
%
% Returns Tp with .S (the template surface) .dk (DK label per vertex) .mri (SCS/NCS or []) .file and, per
% hemisphere h, Tp.(h) = struct(gv, u (unit sphere), x (template coordinates, m), code (depth-dmax tile
% code per vertex), dmax). Cached across calls with the same folder, depth and hemispheres.
%
% See also: rheome.scale.grouptiles, rheome.scale.measure_atlas, rheome.geom.tree
%
% Author: Diellor Basha, 2026

    arguments
        dmax (1,1) double {mustBeInteger, mustBePositive}
        opts.Template char = ''
        opts.Hemis string = ["L" "R"]
    end
    persistent cache key
    d = opts.Template;
    if isempty(d), d = getenv('RHEOME_TEMPLATE'); end
    if isempty(d) && ~isempty(getenv('RHEOME_BRAINSTORM')), d = fullfile(getenv('RHEOME_BRAINSTORM'), 'defaults', 'anatomy', 'ICBM152'); end
    k = string(d) + "|" + dmax + "|" + strjoin(opts.Hemis, ",");
    if isequal(key, k), Tp = cache;  return; end
    f = fullfile(d, 'tess_cortex_pial_low.mat');
    if ~isfile(f)
        error('scale:atlas:template', ['No template cortex at "%s": pass Template=, or set RHEOME_TEMPLATE, or ' ...
              'RHEOME_BRAINSTORM (a Brainstorm checkout with defaults/anatomy/ICBM152).'], f);
    end
    Tp.S = rheome.io.read.surface(f);
    if isempty(Tp.S.Sphere) || isempty(Tp.S.Hemi), error('scale:atlas:template', '%s has no Reg.Sphere or no L/R split', f); end
    A = rheome.io.read.atlas(f);  Tp.dk = i_labels(A(strcmpi({A.Name}, 'Desikan-Killiany')), Tp.S.nV);
    Tp.mri = [];  mf = fullfile(d, 'subjectimage_T1.mat');  if isfile(mf), Tp.mri = load(mf, 'SCS', 'NCS'); end
    Tp.file = f;
    for h = opts.Hemis(:)'
        Sh = rheome.utils.hemisphere(Tp.S, [], char(h));
        Tr = rheome.geom.tree(Sh, MaxDepth=dmax);  G = rheome.geom.tiles(Tr, Sh, dmax);
        if ~G.heap || any(G.depth ~= dmax), error('scale:atlas:heap', 'template hemisphere %s: tree not full to depth %d', h, dmax); end
        gv = Sh.GlobalVertices(:);  u = double(Tp.S.Sphere(gv, :));
        Tp.(h) = struct('gv', gv, 'u', u ./ vecnorm(u, 2, 2), 'x', double(Tp.S.Vertices(gv, :)), ...
                        'code', G.node_id(G.tileOf), 'dmax', dmax);
    end
    cache = Tp;  key = k;
end

function lab = i_labels(A, nV)
% one label per vertex; unknown / medial wall / corpus callosum left empty (as rheome.scale.measure_atlas)
    keep = ~contains(lower(string(A.Label(:))), ["unknown" "medial" "corpus"]);
    L = string(A.Label(:));  L = L(keep);  [mx, k] = max(double(A.Membership(keep, :)), [], 1);
    lab = strings(nV, 1);  lab(mx > 0) = L(k(mx > 0));
end

% Author: Diellor Basha, 2026
