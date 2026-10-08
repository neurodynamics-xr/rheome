function R = rows(db, relation, opts)
% SELECT.ROWS  A relation of the schema, materialised from the store as a MATLAB table.
%
%   R = rheome.select.rows(db, "recording" | "channel" | "band")
%   R = rheome.select.rows(db, "tile" | "tile_cone" | "feature" | "feature_band", Level=L)
%   R = rheome.select.rows(db, "sensor_group")
%   R = rheome.select.rows(db, "feature_group" | "feature_group_band", Level=L)     % internal nodes
%   R = rheome.select.rows(db, "feature_cortex" | "feature_cortex_space", Level=L)  % <store>__cortex.mat
%   R = rheome.select.rows(db, "cortex_node")                                       % likewise, if written
%
% The reference executor's view: exactly the rows a database would hold, in key order.
% feature_band holds only the bands carried at the level (the diagonal).
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        relation (1,1) string
        opts.Level double = []
    end
    meta = db.meta;  g = db.grid;  b = db.bands;  rid = string(db.recording_id);
    L = opts.Level;
    needsLevel = ["tile","tile_cone","feature","feature_band","feature_group","feature_group_band","feature_space","feature_group_space","feature_cortex","feature_cortex_space"];
    if any(relation == needsLevel) && isempty(L)
        error('select:rows:level', 'Relation %s needs Level.', relation);
    end
    if any(relation == ["feature","feature_band","feature_space"]) && L < db.channelLevel
        error('select:level', 'Channel rows are carried from level %d (%g s tiles); level %d is below it. Use Scope="group" or a coarser level.', db.channelLevel, g.tExtent(db.channelLevel+1), L);
    end
    switch relation
        case "recording"
            R = table(rid, string(meta.source), meta.fs, meta.nT, meta.duration, string(meta.units), ...
                      string(i_bank(meta)), string(meta.hash), string(meta.release), string(meta.commit), string(meta.created), ...
                      'VariableNames', {'recording_id','source','fs','n_samples','duration','units','bank','cfg_hash','release','commit','created'});
        case "channel"
            C = meta.C;
            R = table(repmat(rid, C, 1), (1:C)', string(meta.ChannelName(:)), string(meta.ChannelType(:)), meta.iChannel(:), ...
                      'VariableNames', {'recording_id','channel_id','name','type','store_row'});
        case "band"
            nB = height(b);
            if ismember('kind', b.Properties.VariableNames), kind = string(b.kind); else, kind = repmat("octave", nB, 1); end
            R = table(repmat(rid, nB, 1), b.j, kind, b.fLo, b.fHi, b.fCenter, b.fExtent, b.tSupport, b.naturalLevel, b.partial, ...
                      'VariableNames', {'recording_id','band_id','kind','f_lo','f_hi','f_center','f_extent','t_support','natural_level','partial'});
        case "tile"
            K = g.K(L+1);  n = double(rheome.select.level(db, L, 'n'));
            R = table(repmat(rid, K, 1), repmat(L, K, 1), (1:K)', g.tCenter{L+1}(:), repmat(g.tExtent(L+1), K, 1), n(:), ...
                      'VariableNames', {'recording_id','level','k','t_center','t_extent','n'});
        case "tile_cone"
            K = g.K(L+1);  present = g.bandsAt{L+1};  nc = double(rheome.select.level(db, L, 'nCoi'));
            nc = reshape(nc, K, numel(present));
            [kk, bb] = ndgrid(1:K, present);
            R = table(repmat(rid, numel(kk), 1), repmat(L, numel(kk), 1), kk(:), bb(:), nc(:), ...
                      'VariableNames', {'recording_id','level','k','band_id','n_coi'});
            R = sortrows(R, {'k','band_id'});
        case "feature"
            K = g.K(L+1);  C = meta.C;
            sx = rheome.select.level(db, L, 'sumX');  sx2 = rheome.select.level(db, L, 'sumX2');
            am = rheome.select.level(db, L, 'absMax');  mn = rheome.select.level(db, L, 'min');  mx = rheome.select.level(db, L, 'max');
            [kk, cc] = ndgrid(1:K, 1:C);
            R = table(repmat(rid, K*C, 1), cc(:), repmat(L, K*C, 1), kk(:), sx(:), sx2(:), am(:), mn(:), mx(:), ...
                      'VariableNames', {'recording_id','channel_id','level','k','sum_x','sum_x2','abs_max','min','max'});
            R = sortrows(R, {'channel_id','k'});
        case "feature_band"
            K = g.K(L+1);  C = meta.C;  present = g.bandsAt{L+1};  nb = numel(present);
            e = reshape(rheome.select.level(db, L, 'energy'), K, C, nb);  v = reshape(rheome.select.level(db, L, 'envMax'), K, C, nb);
            [kk, cc, bb] = ndgrid(1:K, 1:C, present);
            R = table(repmat(rid, K*C*nb, 1), cc(:), repmat(L, K*C*nb, 1), kk(:), bb(:), e(:), v(:), ...
                      'VariableNames', {'recording_id','channel_id','level','k','band_id','energy','env_max'});
            R = sortrows(R, {'channel_id','k','band_id'});
        case "sensor_group"
            tr = db.tree;
            if isempty(tr), R = table(); return; end
            n = height(tr);  cen = tr.centroid;  if iscell(cen), cen = vertcat(cen{:}); end
            R = table(repmat(rid, n, 1), tr.node_id, tr.parent_id, tr.depth, tr.n_sensors, tr.is_leaf, tr.channel_id, tr.diameter, cen(:,1), cen(:,2), cen(:,3), ...
                      'VariableNames', {'recording_id','group_id','parent_id','depth','n_sensors','is_leaf','channel_id','diameter','centroid_x','centroid_y','centroid_z'});
        case "feature_group"
            K = g.K(L+1);  nodes = db.groupNodes;  nG = numel(nodes);
            sx = rheome.select.level(db, L, 'g_sumX');  sx2 = rheome.select.level(db, L, 'g_sumX2');
            am = rheome.select.level(db, L, 'g_absMax');  mn = rheome.select.level(db, L, 'g_min');  mx = rheome.select.level(db, L, 'g_max');
            [kk, gg] = ndgrid(1:K, nodes);
            R = table(repmat(rid, K*nG, 1), gg(:), repmat(L, K*nG, 1), kk(:), sx(:), sx2(:), am(:), mn(:), mx(:), ...
                      'VariableNames', {'recording_id','group_id','level','k','sum_x','sum_x2','abs_max','min','max'});
            R = sortrows(R, {'group_id','k'});
        case "feature_group_band"
            K = g.K(L+1);  nodes = db.groupNodes;  nG = numel(nodes);  present = g.bandsAt{L+1};  nb = numel(present);
            e = reshape(rheome.select.level(db, L, 'g_energy'), K, nG, nb);  v = reshape(rheome.select.level(db, L, 'g_envMax'), K, nG, nb);
            [kk, gg, bb] = ndgrid(1:K, nodes, present);
            R = table(repmat(rid, K*nG*nb, 1), gg(:), repmat(L, K*nG*nb, 1), kk(:), bb(:), e(:), v(:), ...
                      'VariableNames', {'recording_id','group_id','level','k','band_id','energy','env_max'});
            R = sortrows(R, {'group_id','k','band_id'});
        case "space_band"
            sb = db.sbands;
            if isempty(sb), R = table(); return; end
            n = height(sb);
            R = table(repmat(rid, n, 1), sb.sband_id, string(sb.kind), sb.k_lo, sb.k_hi, sb.k_center, sb.wavelength_lo, sb.wavelength_hi, sb.wavelength_center, sb.n_modes, ...
                      'VariableNames', {'recording_id','sband_id','kind','k_lo','k_hi','k_center','wavelength_lo','wavelength_hi','wavelength_center','n_modes'});
        case "feature_space"
            K = g.K(L+1);  C = meta.C;  cs = db.chanSpace;  nS = numel(cs);
            e = reshape(rheome.select.level(db, L, 's_energy'), K, C, nS);  v = reshape(rheome.select.level(db, L, 's_envMax'), K, C, nS);
            [kk, cc, ss] = ndgrid(1:K, 1:C, cs);
            R = table(repmat(rid, K*C*nS, 1), cc(:), repmat(L, K*C*nS, 1), kk(:), ss(:), e(:), v(:), ...
                      'VariableNames', {'recording_id','channel_id','level','k','sband_id','energy','env_max'});
            R = sortrows(R, {'channel_id','k','sband_id'});
        case "feature_group_space"
            K = g.K(L+1);  sl = db.spaceSlots;  nSl = height(sl);
            e = rheome.select.level(db, L, 'gs_energy');  v = rheome.select.level(db, L, 'gs_envMax');
            kk = repmat((1:K)', nSl, 1);  gg = repelem(sl.group_id, K);  ss = repelem(sl.sband_id, K);
            R = table(repmat(rid, K*nSl, 1), gg, repmat(L, K*nSl, 1), kk, ss, e(:), v(:), ...
                      'VariableNames', {'recording_id','group_id','level','k','sband_id','energy','env_max'});
            R = sortrows(R, {'group_id','k','sband_id'});
        case "cortex_node"
            % ⚠⚠ THE STORE DOES NOT CARRY THIS DIMENSION: cortex_node is rheome.geom.tree's nodes, which live
            % with the SUBJECT'S SURFACE (rheome.load.bases + rheome.geom.tree), not with the recording. The rows come
            % from the cortical sidecar when rheome.flow.cortexfeatures has written one (it records the nodes it
            % summed over); otherwise the empty typed relation -- the columns are real, the rows are
            % elsewhere. ⭐ Adding the relation to rheome.select.schema without this branch is what made
            % tSelectSchema/everyRelationMaterialisesWithItsColumns fail, caught by the full suite
            % rather than by the targeted tests run at the time.
            F = i_cortex(db, "nodes");
            if isempty(F), R = i_empty('cortex_node'); else, R = F.nodes; end
        case "feature_cortex"
            F = i_cortex(db, "lev");
            if isempty(F), R = i_empty('feature_cortex'); return; end
            X = F.lev{L+1};  [nN, K, ~] = size(X);  [nn, kk] = ndgrid(F.ids, 1:K);
            R = table(repmat(rid, nN*K, 1), nn(:), repmat(L, nN*K, 1), kk(:), repmat(F.band_id, nN*K, 1), ...
                      reshape(X(:,:,1), [], 1), reshape(X(:,:,2), [], 1), reshape(X(:,:,3), [], 1), ...
                      'VariableNames', {'recording_id','node_id','level','k','band_id','area_s','on_area_s','energy'});
            R = sortrows(R, {'node_id','k'});
        case "feature_cortex_space"
            F = i_cortex(db, "lev");
            if isempty(F), R = i_empty('feature_cortex_space'); return; end
            X = F.lev{L+1}(:,:,4:end);  [nN, K, nC] = size(X);  [nn, kk, cc] = ndgrid(F.ids, 1:K, 1:nC);
            R = table(repmat(rid, numel(nn), 1), nn(:), repmat(L, numel(nn), 1), kk(:), repmat(F.band_id, numel(nn), 1), ...
                      cc(:), X(:), 'VariableNames', {'recording_id','node_id','level','k','band_id','cband_id','energy'});
            R = sortrows(R, {'node_id','k','cband_id'});
        otherwise
            if any(relation == ["label","measure","txn"])
                error('select:rows:note', ...
                    ['%s is a note relation: it lives in the store''s sidecar file, not in the ' ...
                     'store, and nothing in it merges. Read it with select.%s.'], relation, ...
                     i_reader(relation));
            end
            error('select:rows:relation', 'Unknown relation %s.', relation);
    end
end

function F = i_cortex(db, var)
% The cortical sidecar's variables, or [] when rheome.flow.cortexfeatures has not written one.
    F = [];
    if ~isfield(db, 'file'), return; end
    f = rheome.select.cortexfile(db);
    if ~isfile(f), return; end
    F = load(f, var, 'ids', 'band_id');
end

function R = i_empty(name)
    sc = rheome.select.schema();  f = sc(strcmp({sc.name}, name));
    R = table('Size', [0 numel(f.columns)], 'VariableTypes', repmat({'double'}, 1, numel(f.columns)), ...
              'VariableNames', f.columns);
    for c = intersect(string(f.columns(strcmp(f.types, 'text'))), string(f.columns))
        R.(c) = strings(0, 1);
    end
end

function r = i_reader(relation)
    switch relation
        case "label", r = "labels";  case "measure", r = "measures";  otherwise, r = "labels";
    end
end

function b = i_bank(meta)
% stores built before cfg.Bank existed are Morse stores
    if isfield(meta.cfg, 'Bank'), b = meta.cfg.Bank; else, b = 'morse'; end
end
% Author: Diellor Basha, 2026
