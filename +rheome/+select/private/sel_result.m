function R = sel_result(db, q, L, hits, V, chans, bandsAt)
% SEL_RESULT  Rows for the tiles where hits is true, in key order; duration runs applied.
%   hits, V: [K x nChan x nBand] (nBand = 1 for all-pass stats)
%
% Author: Diellor Basha, 2026

    g = db.grid;  K = g.K(L+1);  rid = string(db.recording_id);
    if isempty(bandsAt), bid = 0; else, bid = bandsAt; end
    [kk, cc, bb] = ndgrid(1:K, chans, bid);
    idx = find(hits(:));                                       % (:) -- a 1 x n x m array would give a row
    tcs = g.tCenter{L+1};
    R = table(repmat(rid, numel(idx), 1), cc(idx), repmat(L, numel(idx), 1), kk(idx), ...
              reshape(tcs(kk(idx)), [], 1), repmat(g.tExtent(L+1), numel(idx), 1), bb(idx), ...
              repmat(string(q.stat), numel(idx), 1), V(idx), ...
              'VariableNames', {'recording_id','channel_id','level','k','t_center','t_extent','band_id','stat','value'});
    if isfield(q, 'scope') && strcmp(q.scope, 'group')
        R.Properties.VariableNames{'channel_id'} = 'group_id';
        R = sortrows(R, {'group_id','band_id','k'});
    else
        R = sortrows(R, {'channel_id','band_id','k'});
    end
    if q.minDuration > 0
        need = ceil(q.minDuration / g.tExtent(L+1));
        R = i_runs(R, need);
    end
end

function R = i_runs(R, need)
% gaps-and-islands: consecutive k within (channel, band); keep runs of >= need tiles
    if isempty(R)
        R.run_id = zeros(0,1);  R.run_start = zeros(0,1);  R.run_length = zeros(0,1);
        return
    end
    unit = R.(R.Properties.VariableNames{2});                   % channel_id or group_id
    grp = unit * 1e9 + R.band_id * 1e6;                         % (unit, band) key
    isNew = [true; diff(grp) ~= 0 | diff(R.k) ~= 1];
    run = cumsum(isNew);
    len = accumarray(run, 1);
    first = accumarray(run, R.k, [], @min);
    keep = len(run) >= need;
    R = R(keep, :);
    run = run(keep);
    [~, ~, rid] = unique(run, 'stable');
    R.run_id = rid;
    R.run_start = first(run);
    R.run_length = len(run);
end
% Author: Diellor Basha, 2026
