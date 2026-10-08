function [V, chans, bandsAt, K] = sel_values(db, L, stat, chans, bandSel, scope, groupSel, depth)
% SEL_VALUES  A statistic at a level as a [K x nChan x nBand] array over the query's scope.
%   bandSel: band ids ([] = all carried at L); returns bandsAt (the ids in the third axis).
%   Derived stats use n from the tile relation. Counts examined tiles in db.cost.
%
% Author: Diellor Basha, 2026

    g = db.grid;  K = g.K(L+1);  C = db.meta.C;
    if nargin < 6 || isempty(scope), scope = 'channel'; end
    if nargin < 7, groupSel = []; end
    if nargin < 8, depth = []; end
    pre = '';
    % ⭐ min and max may reach below the channel floor: they merge by extremum, so the pair
    % on a tile contains every sample in it, which is the min/max mipmap an envelope is drawn
    % from. Every other channel column still starts at the floor.
    % ⚠ AND THE SPECTRAL PEAKS REACH BELOW IT TOO, for a different reason: a peak does not
    % merge, so it is stored at the band's HOME level -- which for alpha is 2 s tiles, far
    % under a 16 s channel floor. The floor is about what the pyramid carries UPWARD; these
    % two families are written where they are computed and are read there.
    if ~strcmp(scope, 'group') && L < db.channelLevel
        envOK  = any(strcmp(stat, {'min','max'})) && L >= db.envelopeLevel;
        peakOK = any(strcmp(stat, {'peakFreq','peakAmp','fftPower'}));
        if ~envOK && ~peakOK
            error('select:level', 'Channel rows are carried from level %d (%g s tiles); level %d is below it. Use Scope="group", a coarser level, or the min/max envelope, which reaches level %d.', ...
                  db.channelLevel, g.tExtent(db.channelLevel+1), L, db.envelopeLevel);
        end
    end
    if strcmp(scope, 'group')
        % the unit axis is the sensor tree's internal nodes; 'chans' become node ids
        nodes = db.groupNodes;
        if isempty(nodes), error('select:groups', 'This store has no sensor groups (built without positions).'); end
        if ~isempty(groupSel), sel = groupSel; else, sel = nodes; end
        if ~isempty(depth), sel = sel(ismember(sel, db.tree.node_id(db.tree.depth == depth))); end
        miss = setdiff(sel, nodes);
        if ~isempty(miss), error('select:groups', 'Node %d is not an internal node with rows.', miss(1)); end
        chans = sel;
        C = numel(nodes);
        col = arrayfun(@(n) find(nodes == n, 1), sel);                      % columns in the group arrays
        pre = 'g_';
    else
        if isempty(chans), chans = 1:C; end
        col = chans;
    end
    present = g.bandsAt{L+1};
    switch stat
        case {'spaceEnergy','spaceEnvMax'}
            nS = height(db.sbands);
            if nS == 0, error('select:space', 'This store has no wavelength axis (built without positions or Space=false).'); end
            if isempty(bandSel), bandsAt = 1:nS; else, bandsAt = bandSel; end
            if any(bandsAt < 1 | bandsAt > nS), error('select:space', 'Spatial band out of range 1..%d.', nS); end
            if strcmp(stat, 'spaceEnergy'), nm = 's_energy'; else, nm = 's_envMax'; end
            if strcmp(scope, 'group')
                % packed slots: (band, node) pairs the diagonal stored
                sl = db.spaceSlots;
                A = rheome.select.level(db, L, ['g' nm]);
                if isempty(groupSel) && ~isempty(depth)
                    % [] means "the nodes at this depth that carry the band"
                    sel = sel(arrayfun(@(n) all(ismember(bandsAt, sl.sband_id(sl.group_id == n))), sel));
                    chans = sel;
                end
                V = zeros(K, numel(sel), numel(bandsAt));
                for i = 1:numel(sel)
                    for j = 1:numel(bandsAt)
                        r = find(sl.group_id == sel(i) & sl.sband_id == bandsAt(j), 1);
                        if isempty(r)
                            error('select:space', 'Node %d (diameter %.0f mm) does not carry spatial band %d (wavelengths from %.0f mm): the position diagonal.', ...
                                  sel(i), 1e3*db.tree.diameter(sel(i)), bandsAt(j), 1e3*db.sbands.wavelength_lo(bandsAt(j)));
                        end
                        V(:, i, j) = A(:, sl.slot(r));
                    end
                end
            else
                cs = db.chanSpace;
                miss = setdiff(bandsAt, cs);
                if ~isempty(miss)
                    error('select:space', 'A single sensor carries only spatial band(s) %s (the position diagonal); band %d is not stored per channel.', mat2str(cs), miss(1));
                end
                A = reshape(rheome.select.level(db, L, nm), K, C, numel(cs));
                slot = arrayfun(@(b) find(cs == b, 1), bandsAt);
                V = double(A(:, col, slot));
            end
        case {'energy','envMax','bandPower'}
            if isempty(bandSel), bandsAt = present; else, bandsAt = bandSel; end
            miss = setdiff(bandsAt, present);
            if ~isempty(miss)
                nl = db.bands.naturalLevel(miss(1));
                error('select:level', 'Band %d is carried from level %d (its natural level); level %d is below it.', miss(1), nl, L);
            end
            slot = arrayfun(@(b) find(present == b, 1), bandsAt);
            if strcmp(stat, 'bandPower')
                A = reshape(rheome.select.level(db, L, [pre 'energy']), K, C, numel(present));
                n = double(rheome.select.level(db, L, 'n'));
                V = A(:, col, slot) ./ max(n, 1);
            else
                A = reshape(rheome.select.level(db, L, [pre stat]), K, C, numel(present));
                V = double(A(:, col, slot));
            end
        case 'nCoi'
            if isempty(bandSel), bandsAt = present; else, bandsAt = bandSel; end
            slot = arrayfun(@(b) find(present == b, 1), bandsAt);
            A = reshape(double(rheome.select.level(db, L, 'nCoi')), K, numel(present));
            V = repmat(reshape(A(:, slot), K, 1, []), 1, numel(col), 1);
        case {'peakFreq','peakAmp','fftPower'}
            % ⚠ ONE LEVEL, AND ONLY ONE. A spectral peak does not merge, so it is stored at
            % the band's HOME level (the first tile at least its support long) and nowhere
            % else. Asking elsewhere is refused with the level that has it.
            if isempty(db.homeAt) || numel(db.homeAt) < L+1
                error('select:peaks', 'This store has no spectral peaks (built with Peaks=false).');
            end
            home = db.homeAt{L+1};
            if isempty(bandSel), bandsAt = home; else, bandsAt = bandSel; end
            miss = setdiff(bandsAt, home);
            if ~isempty(miss)
                hl = db.bands.naturalLevel(miss(1));
                error('select:peaks:level', ...
                      ['Band %d''s spectral peak lives at level %d (its home), not at level %d. ' ...
                       'A peak is not mergeable, so it is stored once.'], miss(1), hl, L);
            end
            if strcmp(scope, 'group')
                error('select:peaks:scope', 'Spectral peaks are per channel; a node is a sum of sensors.');
            end
            nm = struct('peakFreq','pkFreq','peakAmp','pkAmp','fftPower','pkPower');
            A = reshape(rheome.select.level(db, L, nm.(stat)), K, C, numel(home));
            slot = arrayfun(@(b) find(home == b, 1), bandsAt);
            V = double(A(:, col, slot));
        case {'sumX2','absMax','min','max'}
            bandsAt = [];
            V = double(rheome.select.level(db, L, [pre stat]));  V = V(:, col);
        case 'mean'
            bandsAt = [];
            n = double(rheome.select.level(db, L, 'n'));  V = rheome.select.level(db, L, [pre 'sumX']);  V = V(:, col) ./ max(n, 1);
        case 'rms'
            bandsAt = [];
            n = double(rheome.select.level(db, L, 'n'));  V = rheome.select.level(db, L, [pre 'sumX2']);  V = sqrt(V(:, col) ./ max(n, 1));
        case {'sumAbs','sumSqrt','sumX3','sumX4','sumX'}
            bandsAt = [];
            V = double(rheome.select.level(db, L, [pre stat]));  V = V(:, col);
    end
    db.cost('tiles') = db.cost('tiles') + numel(V);
end
% Author: Diellor Basha, 2026
