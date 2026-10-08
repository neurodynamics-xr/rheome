function q = query(opts)
% SELECT.QUERY  A query as a value: stat, operator, threshold, level, scope. Executors interpret it.
%
%   q = rheome.select.query(Stat="envMax", Op=">", Threshold="p95", Level=3, Band=9, Channels=[], MinDuration=0)
%
% Stats: energy, envMax (feature_band, need Band); sumX2, absMax, min, max (feature);
% spaceEnergy, spaceEnvMax (feature_space: Band names a SPATIAL band, sband_id);
% derived: mean = sum_x/n, rms = sqrt(sum_x2/n), bandPower = energy/n.
% Threshold: a number, "pNN" (percentile of the stat over the query's scope at its
% level, resolved per recording), or "a*median".
% Band: one band_id, or [] for every band carried at the level (band stats only).
% Channels: channel_ids, or [] for all. MinDuration in seconds (0 = none).
% Scope: "channel" (rows per channel) | "group" (rows per sensor-tree node, the position
% pyramid); Groups: node ids ([] = all internal nodes), Depth: only nodes at that tree depth.
%
% The struct carries the relation and the pruning rule (design §3.2):
%   .bound  the parent test: 'ge' (parent >= theta), 'le' (parent <= theta), 'none',
%           or a proxy stat with its rule, e.g. rms -> absMax 'ge'
%
% Author: Diellor Basha, 2026

    arguments
        opts.Stat        (1,1) string {mustBeMember(opts.Stat, ["energy","envMax","sumX2","absMax","min","max","mean","rms","bandPower","nCoi","spaceEnergy","spaceEnvMax"])}
        opts.Op          (1,1) string {mustBeMember(opts.Op, [">",">=","<","<="])} = ">"
        opts.Threshold   = "p95"
        opts.Level       (1,1) double {mustBeInteger, mustBeNonnegative}
        opts.Band        double = []
        opts.Channels    double = []
        opts.MinDuration (1,1) double {mustBeNonnegative} = 0
        opts.Scope       (1,1) string {mustBeMember(opts.Scope, ["channel","group"])} = "channel"
        opts.Groups      double = []
        opts.Depth       double = []
    end
    q = struct();
    q.stat = char(opts.Stat);  q.op = char(opts.Op);  q.threshold = opts.Threshold;
    q.level = opts.Level;  q.band = opts.Band(:)';  q.channels = opts.Channels(:)';  q.minDuration = opts.MinDuration;
    q.scope = char(opts.Scope);  q.groups = opts.Groups(:)';  q.depth = opts.Depth;
    if strcmp(q.scope, 'group') && ~isempty(q.channels)
        error('select:query:scope', 'Channels selects channel rows; use Groups or Depth for group scope.');
    end
    switch q.stat
        case {'energy','envMax'},        q.relation = 'feature_band';  q.column = i_col(q.stat);
        case {'sumX2','absMax','min','max'}, q.relation = 'feature';  q.column = i_col(q.stat);
        case 'nCoi',                     q.relation = 'tile_cone';     q.column = 'n_coi';
        case {'mean','rms'},             q.relation = 'feature';       q.column = q.stat;   % derived
        case 'bandPower',                q.relation = 'feature_band';  q.column = q.stat;   % derived
        case {'spaceEnergy','spaceEnvMax'}, q.relation = 'feature_space';  q.column = i_col(q.stat);
    end
    if strcmp(q.scope, 'group') && ~strcmp(q.relation, 'tile_cone')
        q.relation = strrep(q.relation, 'feature', 'feature_group');      % feature -> feature_group, feature_band -> feature_group_band
    end
    q.derived = any(strcmp(q.stat, {'mean','rms','bandPower'}));
    q.bandStat = any(strcmp(q.relation, {'feature_band','tile_cone','feature_group_band','feature_space','feature_group_space'}));
    q.spaceStat = any(strcmp(q.stat, {'spaceEnergy','spaceEnvMax'}));
    greater = any(strcmp(q.op, {'>','>='}));
    % pruning rule: which parent statistic bounds the child, and how
    switch q.stat
        case {'energy','sumX2','nCoi','envMax','absMax','max','spaceEnergy','spaceEnvMax'}
            if greater, q.bound = struct('stat', q.stat, 'rule', 'ge'); else, q.bound = struct('stat', '', 'rule', 'none'); end
        case 'min'
            if ~greater, q.bound = struct('stat', 'min', 'rule', 'le'); else, q.bound = struct('stat', '', 'rule', 'none'); end
        case {'mean','rms'}
            if greater, q.bound = struct('stat', 'absMax', 'rule', 'ge'); else, q.bound = struct('stat', '', 'rule', 'none'); end
        case 'bandPower'
            if greater, q.bound = struct('stat', 'envMax', 'rule', 'ge_sqrt'); else, q.bound = struct('stat', '', 'rule', 'none'); end
    end
end

function c = i_col(stat)
    switch stat
        case 'energy', c = 'energy';  case 'envMax', c = 'env_max';  case 'sumX2', c = 'sum_x2';
        case 'absMax', c = 'abs_max'; case 'min', c = 'min';         case 'max', c = 'max';
        case 'spaceEnergy', c = 'energy';  case 'spaceEnvMax', c = 'env_max';
    end
end
% Author: Diellor Basha, 2026
