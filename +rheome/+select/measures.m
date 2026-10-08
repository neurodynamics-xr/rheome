function [M, T] = measures(db, opts)
% SELECT.MEASURES  Read measurements (and their transactions), filtered.
%
%   M = rheome.select.measures(db)
%   [M, T] = rheome.select.measures(db, Kind="vortex_count", Level=3, Scope="group", Units=[2 3])
%   M = rheome.select.measures(db, Kind="vortex_count", Min=2, Window=[120 180])
%
% Window filters on the node's own span, so a measurement on a coarse node is returned for
% any window it overlaps -- the same rule rheome.select.derive uses, for the same reason.
%
% See also: rheome.select.measure, rheome.select.context
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        opts.Kind   string = string.empty
        opts.Level  double = []
        opts.Scope  string = string.empty
        opts.Units  double = []
        opts.Band   double = []
        opts.Txn    double = []
        opts.Min    double = []
        opts.Max    double = []
        opts.Window double = []
    end
    [~, T, M] = sel_labelio(db);
    keep = true(height(M), 1);
    if ~isempty(opts.Kind),  keep = keep & ismember(M.kind, opts.Kind); end
    if ~isempty(opts.Level), keep = keep & ismember(M.level, opts.Level); end
    if ~isempty(opts.Scope), keep = keep & ismember(M.scope, opts.Scope); end
    if ~isempty(opts.Units), keep = keep & ismember(M.unit_id, opts.Units); end
    if ~isempty(opts.Band),  keep = keep & ismember(M.band_id, opts.Band); end
    if ~isempty(opts.Txn),   keep = keep & ismember(M.txn_id, opts.Txn); end
    if ~isempty(opts.Min),   keep = keep & M.value >= opts.Min; end
    if ~isempty(opts.Max),   keep = keep & M.value <= opts.Max; end
    if ~isempty(opts.Window)
        g = db.grid;  w = sort(double(opts.Window(:)'));
        ext = g.tExtent(M.level + 1)';
        lo = (M.k - 1) .* ext;  hi = min(lo + ext, db.meta.duration);
        keep = keep & hi > w(1) & lo < w(2);
    end
    M = M(keep, :);
    if nargout > 1, T = T(ismember(T.txn_id, unique(M.txn_id)), :); end
end
% Author: Diellor Basha, 2026
