function txt = ddl()
% SELECT.DDL  The schema as Postgres DDL, with the indexes the query plans assume.
%
%   txt = rheome.select.ddl()
%
% Author: Diellor Basha, 2026

    S = rheome.select.schema();
    lines = {};
    for r = S
        cols = cellfun(@(c, t) sprintf('  %s %s', c, t), r.columns, r.types, 'UniformOutput', false);
        cols{end+1} = sprintf('  PRIMARY KEY (%s)', strjoin(r.key, ', ')); %#ok<AGROW>
        fk = fieldnames(r.foreign);
        for i = 1:numel(fk)
            switch fk{i}
                case 'level_k', cols{end+1} = '  FOREIGN KEY (recording_id, level, k) REFERENCES tile (recording_id, level, k)'; %#ok<AGROW>
                case 'channel_id', cols{end+1} = '  FOREIGN KEY (recording_id, channel_id) REFERENCES channel (recording_id, channel_id)'; %#ok<AGROW>
                case 'band_id', cols{end+1} = '  FOREIGN KEY (recording_id, band_id) REFERENCES band (recording_id, band_id)'; %#ok<AGROW>
                case 'group_id', cols{end+1} = '  FOREIGN KEY (recording_id, group_id) REFERENCES sensor_group (recording_id, group_id)'; %#ok<AGROW>
                case 'sband_id', cols{end+1} = '  FOREIGN KEY (recording_id, sband_id) REFERENCES space_band (recording_id, sband_id)'; %#ok<AGROW>
                otherwise, cols{end+1} = sprintf('  FOREIGN KEY (%s) REFERENCES %s (%s)', fk{i}, r.foreign.(fk{i}), fk{i}); %#ok<AGROW>
            end
        end
        lines{end+1} = sprintf('CREATE TABLE %s (\n%s\n);', r.name, strjoin(cols, sprintf(',\n'))); %#ok<AGROW>
    end
    lines{end+1} = 'CREATE INDEX feature_band_env ON feature_band (recording_id, band_id, level, env_max);';
    lines{end+1} = 'CREATE INDEX feature_band_energy ON feature_band (recording_id, band_id, level, energy);';
    lines{end+1} = 'CREATE INDEX feature_absmax ON feature (recording_id, level, abs_max);';
    lines{end+1} = 'CREATE INDEX feature_group_band_env ON feature_group_band (recording_id, band_id, level, env_max);';
    lines{end+1} = 'CREATE INDEX sensor_group_parent ON sensor_group (recording_id, parent_id);';
    lines{end+1} = 'CREATE INDEX feature_space_energy ON feature_space (recording_id, sband_id, level, energy);';
    lines{end+1} = 'CREATE INDEX label_tile ON label (recording_id, level, k, kind);';
    lines{end+1} = 'CREATE INDEX measure_node ON measure (recording_id, kind, level, k);';
    lines{end+1} = 'CREATE INDEX measure_unit ON measure (recording_id, scope, unit_id, level, k);';
    lines{end+1} = 'CREATE INDEX measure_value ON measure (recording_id, kind, value);';
    lines{end+1} = 'CREATE UNIQUE INDEX txn_content ON txn (recording_id, content_hash);';
    txt = strjoin(lines, sprintf('\n\n'));
end
% Author: Diellor Basha, 2026
