function [L, T] = labels(db, opts)
% SELECT.LABELS  Read labels (and transactions) of a store, filtered.
%
%   [L, T] = rheome.select.labels(db)
%   L = rheome.select.labels(db, Kind="bad", Level=3, Channels=[1 2], Txn=id)
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        opts.Kind     string = string.empty
        opts.Level    double = []
        opts.Channels double = []
        opts.Txn      double = []
    end
    [L, T] = sel_labelio(db);
    keep = true(height(L), 1);
    if ~isempty(opts.Kind),     keep = keep & ismember(L.kind, opts.Kind); end
    if ~isempty(opts.Level),    keep = keep & ismember(L.level, opts.Level); end
    if ~isempty(opts.Channels), keep = keep & ismember(L.channel_id, opts.Channels); end
    if ~isempty(opts.Txn),      keep = keep & ismember(L.txn_id, opts.Txn); end
    L = L(keep, :);
end
% Author: Diellor Basha, 2026
