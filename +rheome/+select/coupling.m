function T = coupling(db, opts)
% SELECT.COUPLING  Phase-amplitude coupling per tile, derived from the stored sums.
%
%   T = rheome.select.coupling(db)                                  % every pair, at its own level
%   T = rheome.select.coupling(db, Level=8, Channels=66, Window=[200 320])
%   T = rheome.select.coupling(db, Pairs=[7 4])                     % slow band 7, fast band 4
%
% One row per (channel, tile, pair): the mean vector length, the preferred phase, and the
% amplitude sum behind them.
%
% ⭐ ASK AT ANY LEVEL AT OR ABOVE THE PAIR'S OWN. The store keeps the sums where the pair was
% accumulated -- the slow band's level, 22.6 of its cycles wide. Both are sums, so a coarser
% estimate is the sum of the children's sums: exact, not an average of ratios, and with more
% cycles behind it. Below that level there is nothing to sum, and the refusal says so.
%
% ⚠ MVL IS NOT COMPARABLE ACROSS SAMPLE COUNTS BY ITSELF. |sum|/sum is bounded by 1 and
% biased upward by short windows; on this tiling every pair at its own level has the same
% 22.6 slow cycles, which is what makes the bias a constant of the bank rather than a
% function of frequency. Say which level a number came from.
%
% COLUMNS
%   recording_id, channel_id, level, k, t_center, t_extent
%   slow_band, fast_band, slow_lo, slow_hi, fast_lo, fast_hi, octaves
%   mvl          |sum A exp(i phi)| / sum A      (Canolty's mean vector length)
%   pref_phase   angle of that sum, radians, the slow ANALYTIC phase where fast amplitude peaks
%   amp          sum A over the tile (the normaliser, kept so rows can be re-merged)
%   n_slow_cycles  how many slow cycles the estimate rests on
%
% See also: rheome.ingest.couple, rheome.select.derive, rheome.select.ladder
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        opts.Level    double = []
        opts.Channels double = []
        opts.Pairs    double = []
        opts.Window   double = []
    end
    if ~isfield(db.grid, 'pairsAt') || isempty(db.grid.pairsAt)
        error('select:coupling:none', ...
              'This store has no coupling sums (built with Couple=false).');
    end
    g = db.grid;  b = db.bands;  rid = string(db.recording_id);
    chans = opts.Channels;  if isempty(chans), chans = 1:db.meta.C; end
    T = table();
    for home = 0:g.Lmax
        pr = g.pairsAt{home+1};
        if isempty(pr), continue; end
        keep = true(size(pr, 1), 1);
        if ~isempty(opts.Pairs)
            keep = ismember(pr, opts.Pairs, 'rows');
            if ~any(keep), continue; end
        end
        L = home;  if ~isempty(opts.Level), L = opts.Level; end
        if L < home
            error('select:coupling:level', ...
                  ['Pair (%d,%d) is accumulated at level %d, the slow band''s own; level %d ' ...
                   'is below it and there is nothing to sum.'], pr(1,1), pr(1,2), home, L);
        end
        % ⚠ the store's shape is [K x C x nP], like every other band-resolved array
        V = rheome.select.level(db, home, 'cpVec');
        A = rheome.select.level(db, home, 'cpAmp');
        V = V(:, chans, keep);  A = A(:, chans, keep);
        p = pr(keep, :);
        r = 2^(L - home);
        if r > 1                                              % ⭐ the merge: sums of sums
            K2 = g.K(L+1);
            V = i_fold(V, r, K2);  A = i_fold(A, r, K2);
        end
        K = size(V, 1);
        ext = g.tExtent(L+1);  tc = g.tCenter{L+1}(:);
        lo = (0:K-1)' * ext;  hi = min(lo + ext, db.meta.duration);
        sel = true(K, 1);
        if ~isempty(opts.Window)
            w = sort(double(opts.Window(:)'));
            sel = hi > w(1) & lo < w(2);
        end
        [kk, cc, jj] = ndgrid(find(sel), 1:numel(chans), 1:size(p, 1));
        idx = sub2ind(size(V), kk(:), cc(:), jj(:));
        v = double(V(idx));  a = double(A(idx));
        n = numel(v);
        R = table(repmat(rid, n, 1), reshape(chans(cc(:)), n, 1), repmat(L, n, 1), kk(:), ...
                  tc(kk(:)), repmat(ext, n, 1), p(jj(:), 1), p(jj(:), 2), ...
                  b.fLo(p(jj(:), 1)), b.fHi(p(jj(:), 1)), b.fLo(p(jj(:), 2)), b.fHi(p(jj(:), 2)), ...
                  log2(b.fCenter(p(jj(:), 2)) ./ b.fCenter(p(jj(:), 1))), ...
                  abs(v) ./ max(a, realmin), angle(v), a, ext * b.fCenter(p(jj(:), 1)), ...
                  'VariableNames', {'recording_id','channel_id','level','k','t_center','t_extent', ...
                                    'slow_band','fast_band','slow_lo','slow_hi','fast_lo','fast_hi', ...
                                    'octaves','mvl','pref_phase','amp','n_slow_cycles'});
        T = [T; R];                                            %#ok<AGROW>
    end
    if ~isempty(T), T = sortrows(T, {'channel_id','slow_band','fast_band','k'}); end
end

% Fold r consecutive tiles into one by summing -- the whole point of an additive accumulator.
function Y = i_fold(X, r, K2)
    [K, nC, nP] = size(X);
    pad = r * K2 - K;
    if pad > 0, X(K+1:r*K2, :, :) = 0; end
    Y = reshape(sum(reshape(X(1:r*K2, :, :), r, K2, nC, nP), 1), K2, nC, nP);
end
% Author: Diellor Basha, 2026
