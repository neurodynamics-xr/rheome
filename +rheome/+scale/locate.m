function L = locate(protoDir, sub)
% SCALE.LOCATE  Find one subject's files inside an unzipped Brainstorm protocol.
%
%   L = rheome.scale.locate(protoDir)          % the only subject under data/
%   L = rheome.scale.locate(protoDir, 'sub-01')
%
% The expected protocol layout (one subject per protocol, Brainstorm's own folder names):
%   anat/<sub>/tess_cortex_pial_low.mat                    the cortex the leadfield is built on
%   data/<sub>/@raw<sub>_<ses>_task-rest_run-NN_..._high/  channel, headmodel_surf_os_meg,
%                                                           noisecov_full, data_0raw_*.mat + .bst
%   data/<sub>/@raw<sub>_<ses>_task-noise..._high/         the same-session empty room
%
% ⭐ THE FIRST REST RUN (sorted by name) THAT CARRIES A HEADMODEL is used. Some datasets have two
% rest runs per session; one run keeps the duration comparable with single-run datasets.
% ⚠ The headmodel's SurfaceFile is checked against the cortex: a leadfield on another surface
% would silently mis-index every vertex.
%
% Returns L with .sub .cortex .restDir .restRaw .channel .headmodel .noisecov .noiseRaw
% (.noiseRaw '' if the protocol holds no noise run).
%
% See also: rheome.scale.importsubject, rheome.scale.run
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(sub)
        d = dir(fullfile(protoDir, 'data', 'sub-*'));  d = d([d.isdir]);
        if numel(d) ~= 1, error('scale:locate:subject', 'expected one subject under %s/data, found %d', protoDir, numel(d)); end
        sub = d.name;
    end
    L.sub = char(sub);
    L.cortex = fullfile(protoDir, 'anat', L.sub, 'tess_cortex_pial_low.mat');
    i_need(L.cortex);

    runs = dir(fullfile(protoDir, 'data', L.sub, '@raw*task-rest*'));
    runs = runs([runs.isdir]);  [~, o] = sort({runs.name});  runs = runs(o);
    L.restDir = '';
    for r = runs(:)'
        p = fullfile(r.folder, r.name);
        if ~isempty(dir(fullfile(p, 'headmodel_surf_os_meg*.mat'))), L.restDir = p; break; end
    end
    if isempty(L.restDir), error('scale:locate:rest', 'no rest run with a headmodel under %s', fullfile(protoDir,'data',L.sub)); end
    L.restRaw   = i_one(L.restDir, 'data_0raw_*.mat');
    L.channel   = i_one(L.restDir, 'channel_*.mat');
    L.headmodel = i_one(L.restDir, 'headmodel_surf_os_meg*.mat');
    L.noisecov  = fullfile(L.restDir, 'noisecov_full.mat');  i_need(L.noisecov);

    hm = builtin('load', L.headmodel, 'SurfaceFile');
    [~, hs] = fileparts(char(hm.SurfaceFile));
    if ~strcmp(hs, 'tess_cortex_pial_low')
        error('scale:locate:surface', 'headmodel is on %s, not tess_cortex_pial_low', hm.SurfaceFile);
    end

    nz = dir(fullfile(protoDir, 'data', L.sub, '@raw*task-noise*'));  nz = nz([nz.isdir]);
    L.noiseRaw = '';
    if ~isempty(nz), L.noiseRaw = i_one(fullfile(nz(1).folder, nz(1).name), 'data_0raw_*.mat'); end
end

function f = i_one(d, pat)
    g = dir(fullfile(d, pat));  g = g(~startsWith({g.name}, '._'));
    if isempty(g), error('scale:locate:missing', 'no %s in %s', pat, d); end
    [~, o] = sort([g.datenum], 'descend');  f = fullfile(g(o(1)).folder, g(o(1)).name);
end

function i_need(f)
    if ~exist(f, 'file'), error('scale:locate:missing', 'missing %s', f); end
end

% Author: Diellor Basha, 2026
