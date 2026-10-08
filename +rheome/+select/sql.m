function txt = sql(q, db)
% SELECT.SQL  A query as Postgres SQL: the same semi-join descent the MATLAB executor runs.
%
%   txt = rheome.select.sql(q)          levels from :top down to q.level (parameters :rec, :top)
%   txt = rheome.select.sql(q, db)      with the store's recording id and top level filled in
%
% Renders, in order: a threshold CTE (percentile_cont or a literal), one survivor CTE per
% level from the top down (a semi-join on the parent key (k - 1) / 2 + 1), the final
% SELECT at q.level joined to tile for centre and extent, and for MinDuration a
% gaps-and-islands wrapper. A planner with an index on the column may replace the chain
% by an index scan; the result is the same relation.
%
% Author: Diellor Basha, 2026

    if nargin < 2, rec = ':rec';  top = ':top';  topL = [];  ext = ':t_extent';
    else, rec = sprintf('''%s''', db.recording_id);  topL = db.grid.Lmax;  top = num2str(topL);  ext = num2str(db.grid.tExtent(q.level+1)); end
    rel = q.relation;  col = i_expr(q);
    scope = sprintf('recording_id = %s', rec);
    bandKey = 'band_id';  if q.spaceStat, bandKey = 'sband_id'; end
    if q.bandStat && ~isempty(q.band), scope = [scope sprintf(' AND %s IN (%s)', bandKey, strjoin(string(q.band), ', '))]; end
    if ~isempty(q.channels), scope = [scope sprintf(' AND channel_id IN (%s)', strjoin(string(q.channels), ', '))]; end
    unit = 'channel_id';
    if strcmp(q.scope, 'group')
        unit = 'group_id';
        if ~isempty(q.groups), scope = [scope sprintf(' AND group_id IN (%s)', strjoin(string(q.groups), ', '))]; end
        if ~isempty(q.depth), scope = [scope sprintf(' AND group_id IN (SELECT group_id FROM sensor_group WHERE recording_id = %s AND depth = %d)', rec, q.depth)]; end
    end
    joinTile = q.derived && ~strcmp(q.relation, 'tile_cone');
    from = rel;
    if joinTile, from = sprintf('%s JOIN tile USING (recording_id, level, k)', rel); end
    % threshold
    if isnumeric(q.threshold)
        thr = sprintf('thr AS (SELECT %.15g AS theta)', q.threshold);
    else
        t = char(string(q.threshold));
        tok = regexp(t, '^p(\d+(\.\d+)?)$', 'tokens', 'once');
        if ~isempty(tok)
            thr = sprintf('thr AS (SELECT percentile_cont(%g) WITHIN GROUP (ORDER BY %s) AS theta FROM %s WHERE %s AND level = %d)', ...
                          str2double(tok{1})/100, col, from, scope, q.level);
        else
            tok = regexp(t, '^(\d+(\.\d+)?)\*median$', 'tokens', 'once');
            thr = sprintf('thr AS (SELECT %s * percentile_cont(0.5) WITHIN GROUP (ORDER BY %s) AS theta FROM %s WHERE %s AND level = %d)', ...
                          tok{1}, col, from, scope, q.level);
        end
    end
    ctes = {thr};
    % survivor chain
    keyCols = [unit ', k'];  if q.bandStat, keyCols = [unit ', ' bandKey ', k']; end
    if strcmp(q.relation, 'tile_cone'), keyCols = 'band_id, k'; end
    if ~strcmp(q.bound.rule, 'none')
        [bcol, bexpr] = i_bound(q);
        if ~isempty(topL), levels = topL:-1:q.level+1; else, levels = []; end
        if isempty(levels)
            ctes{end+1} = sprintf('s_top AS (SELECT %s FROM %s WHERE %s AND level = %s AND %s)', keyCols, bcol, scope, top, bexpr);
            ctes{end+1} = sprintf('-- repeat for each level from :top - 1 down to %d: s_L AS (SELECT f.%s FROM %s f JOIN s_(L+1) p ON <parent key> WHERE f.level = L AND %s)', q.level + 1, strrep(keyCols, ', ', ', f.'), bcol, bexpr);
            prev = 's_top';
        else
            prev = '';
            for L = levels
                name = sprintf('s%d', L);
                if isempty(prev)
                    ctes{end+1} = sprintf('%s AS (SELECT %s FROM %s WHERE %s AND level = %d AND %s)', name, keyCols, bcol, scope, L, bexpr); %#ok<AGROW>
                else
                    ctes{end+1} = sprintf('%s AS (SELECT f.%s FROM %s f JOIN %s p ON %s WHERE f.%s AND f.level = %d AND f.%s)', ...
                                          name, strrep(keyCols, ', ', ', f.'), bcol, prev, i_parentjoin(keyCols), strrep(scope, ' AND ', ' AND f.'), L, bexpr); %#ok<AGROW>
                end
                prev = name;
            end
        end
        semi = sprintf(' JOIN %s p ON %s', prev, i_parentjoin(keyCols));
    else
        semi = '';
    end
    bandCol = bandKey;  if ~q.bandStat, bandCol = '0 AS band_id'; end
    chanCol = unit;  if strcmp(q.relation, 'tile_cone'), chanCol = '0 AS channel_id'; end
    hits = sprintf(['hits AS (SELECT f.recording_id, f.%s, f.level, f.k, t.t_center, t.t_extent, f.%s, ''%s'' AS stat, %s AS value ' ...
                    'FROM %s f JOIN tile t ON t.recording_id = f.recording_id AND t.level = f.level AND t.k = f.k%s ' ...
                    'WHERE f.%s AND f.level = %d AND %s %s (SELECT theta FROM thr))'], ...
                   chanCol, bandCol, q.stat, i_qualified(col), rel, semi, strrep(scope, ' AND ', ' AND f.'), q.level, i_qualified(col), q.op);
    ctes{end+1} = hits;
    if q.minDuration > 0
        need = sprintf('CEIL(%g / %s)', q.minDuration, ext);
        ctes{end+1} = sprintf('runs AS (SELECT h.*, h.k - ROW_NUMBER() OVER (PARTITION BY %s, %s ORDER BY k) AS grp FROM hits h)', unit, bandKey);
        % window functions cannot appear in WHERE, so the run columns are computed in a
        % subquery and filtered outside it
        final = sprintf(['SELECT * FROM (SELECT r.*, COUNT(*) OVER (PARTITION BY %s, %s, grp) AS run_length, ' ...
                         'MIN(k) OVER (PARTITION BY %s, %s, grp) AS run_start, DENSE_RANK() OVER (ORDER BY %s, %s, grp) AS run_id ' ...
                         'FROM runs r) x WHERE run_length >= %s ORDER BY %s, %s, k'], unit, bandKey, unit, bandKey, unit, bandKey, need, unit, bandKey);
    else
        final = sprintf('SELECT * FROM hits ORDER BY %s, %s, k', unit, bandKey);
    end
    txt = sprintf('WITH %s\n%s;', strjoin(ctes, sprintf(',\n')), final);
end

function e = i_expr(q)
    switch q.stat
        case 'mean',      e = 'sum_x / NULLIF(n, 0)';
        case 'rms',       e = 'sqrt(sum_x2 / NULLIF(n, 0))';
        case 'bandPower', e = 'energy / NULLIF(n, 0)';
        otherwise,        e = q.column;
    end
end

function e = i_qualified(col)
% qualify bare column names with f. (the derived expressions reference tile columns too)
    e = regexprep(col, '\<(sum_x2|sum_x|abs_max|energy|env_max|n_coi|min|max)\>', 'f.$1');
    e = regexprep(e, '\<n\>(?!_)', 't.n');
end

function [rel, expr] = i_bound(q)
    switch q.bound.stat
        case 'energy',  rel = 'feature_band';  c = 'energy';
        case 'envMax',  rel = 'feature_band';  c = 'env_max';
        case 'sumX2',   rel = 'feature';       c = 'sum_x2';
        case 'absMax',  rel = 'feature';       c = 'abs_max';
        case 'max',     rel = 'feature';       c = 'max';
        case 'min',     rel = 'feature';       c = 'min';
        case 'nCoi',    rel = 'tile_cone';     c = 'n_coi';
        case 'spaceEnergy', rel = 'feature_space';  c = 'energy';
        case 'spaceEnvMax', rel = 'feature_space';  c = 'env_max';
    end
    if strcmp(q.scope, 'group') && ~strcmp(rel, 'tile_cone'), rel = strrep(rel, 'feature', 'feature_group'); end
    switch q.bound.rule
        case 'ge',      expr = sprintf('%s >= (SELECT theta FROM thr)', c);
        case 'le',      expr = sprintf('%s <= (SELECT theta FROM thr)', c);
        case 'ge_sqrt', expr = sprintf('%s * %s >= (SELECT theta FROM thr)', c, c);
    end
end

function j = i_parentjoin(keyCols)
    keys = strsplit(keyCols, ', ');
    parts = cell(1, numel(keys));
    for i = 1:numel(keys)
        if strcmp(keys{i}, 'k'), parts{i} = 'p.k = (f.k - 1) / 2 + 1';
        else, parts{i} = sprintf('p.%s = f.%s', keys{i}, keys{i}); end
    end
    j = strjoin(parts, ' AND ');
end
% Author: Diellor Basha, 2026
