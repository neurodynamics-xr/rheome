function C = catalog(root)
% SELECT.CATALOG  Every tile store under +data (or a folder), one row per store, with what a cohort needs.
%
%   C = rheome.select.catalog()          under rheome.load.root()
%   C = rheome.select.catalog(folder)    under a folder of datasets
%
% Columns: recording_id, dataset, file, store (default | native), bank, voices, frame_floor,
% fs, n_samples, duration, channels, cfg_hash, created. One row per store, so a dataset
% appears once per (store, bank, config). The relational catalogue the spec's SQL table
% stands for; subject metadata joins here when a source for it exists (inventory §1.2:
% none is stored today).
%
% Author: Diellor Basha, 2026

    if nargin < 1, root = rheome.load.root(); end
    d = dir(fullfile(root, '*', 'ingest__*.mat'));
    d = d(~endsWith({d.name}, {'__labels.mat', '__cortex.mat'}));          % the sidecars, not stores
    rows = cell(numel(d), 13);
    for i = 1:numel(d)
        f = fullfile(d(i).folder, d(i).name);
        m = matfile(f);  meta = m.meta;
        [~, ds] = fileparts(d(i).folder);
        store = "default";  if contains(d(i).name, '__native'), store = "native"; end
        rows(i, :) = {string(meta.name), string(ds), string(f), store, string(i_bank(meta)), meta.cfg.VoicesPerOctave, meta.cfg.FrameFloor, meta.fs, meta.nT, meta.duration, meta.C, string(meta.hash), string(meta.created)};
    end
    C = cell2table(rows, 'VariableNames', {'recording_id','dataset','file','store','bank','voices','frame_floor','fs','n_samples','duration','channels','cfg_hash','created'});
    if ~isempty(C), C = sortrows(C, {'dataset','store','bank'}); end
end

function b = i_bank(meta)
% stores built before cfg.Bank existed are Morse stores
    if isfield(meta.cfg, 'Bank'), b = meta.cfg.Bank; else, b = 'morse'; end
end
% Author: Diellor Basha, 2026
