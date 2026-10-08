function [P, pairs] = couple(x, fb, bands, g, opts)
% INGEST.COUPLE  Phase-amplitude coupling accumulators per tile: sums, so they MERGE.
%
%   [P, pairs] = rheome.ingest.couple(x, fb, bands, g)
%   P = rheome.ingest.couple(x, fb, bands, g, Pairs=[7 3], PhaseBins=18)
%
% ⭐ COUPLING IS A SUM, AND THAT CHANGES EVERYTHING. The loudest frequency of a tile does not
% merge; the coupling of two bands does, because every quantity it needs is additive:
%
%     vec(tile) = sum over the tile of  A_fast * exp(i * phi_slow)      (complex)
%     amp(tile) = sum over the tile of  A_fast
%     hist(tile, bin) = the same amplitude sum, split by slow-phase bin
%
% A parent's sums are its children's sums, so the vector length |vec|/amp (Canolty's mean
% vector length), the preferred phase angle(vec) and Tort's modulation index from the
% histogram are all DERIVED at read time, at any level, exactly. Coarser levels give more
% cycles and less variance for nothing.
%
% ⭐ THE SLOW BAND PICKS THE WINDOW. Coupling is estimated over slow cycles, so a pair is
% accumulated at the SLOW band's own level -- where the constant-Q tiling gives 22.6 cycles
% of it, the same number in every band. A fixed window would give 181 cycles at 90 Hz and 1
% at 0.35 Hz, and the estimator's bias would differ by two orders of magnitude across the
% ladder. This is the measurement the tiling was made for.
%
% ⚠ THE FAST BAND PICKS THE RATE. A_fast * exp(i phi_slow) is a real envelope modulated at
% the slow frequency, so the product must be formed at twice the slow band's TOP edge, not
% at the slow band's own rate. Both bands are therefore evaluated on ONE grid -- the fast
% band's -- which on this bank carries every pair two octaves apart or more (measured: 1.2x
% margin at exactly two octaves). Adjacent octaves are refused by name.
%
% ⚠ AND THE PHASE MUST BE RE-MODULATED. The sub-band inverse FFT places a band's bins at the
% start of a short array, which demodulates it by its first bin: the phase that comes back
% rotates against the true one at that frequency, and a sum of phasors built from it is
% destroyed rather than merely offset. Each band is multiplied back by exp(i 2 pi f0 t)
% before any phase is taken.
%
% INPUTS
%   x      [nT x 1] one channel
%   fb     the timefilterbank the store was built with
%   bands  the band table          g  the grid
%   Pairs  [nP x 2] band ids (slow, fast), or "auto" (every octave pair >= 2 octaves apart)
%   PhaseBins  bins of the slow phase for the modulogram (0 = do not accumulate one)
%
% OUTPUT (struct array P), one entry per SLOW level:
%   .level .pairs [nP x 2] .vec [K x nP] complex .amp [K x nP] .hist [K x nP x nBins]
%
% See also: rheome.ingest.peaks, rheome.select.derive, rheome.ingest.bank
%
% Author: Diellor Basha, 2026

    arguments
        x double
        fb (1,1) rheome.timefilterbank
        bands table
        g (1,1) struct
        opts.Pairs = "auto"
        opts.PhaseBins (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        opts.MinOctaves (1,1) double = 2
    end
    pairs = i_pairs(bands, opts.Pairs, opts.MinOctaves);
    P = struct('level', {}, 'pairs', {}, 'vec', {}, 'amp', {}, 'hist', {});
    if isempty(pairs), return; end

    x = double(x(:));  N = g.nT;  fs = g.fs;
    X = fft(x);
    nb = floor(N/2) + 1;
    f = (0:nb-1)' * fs / N;
    B = bins(fb);                                             % [M x 2] per member
    band = i_bandbins(fb, bands, B, f);                       % per band: bins and composite gain

    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for L = unique(bands.naturalLevel(pairs(:,1)))'
        k = find(bands.naturalLevel(pairs(:,1)) == L)';
        K = g.K(L+1);  tile = g.tExtent(L+1);
        e = struct('level', L, 'pairs', pairs(k, :), 'vec', complex(zeros(K, numel(k))), ...
                   'amp', zeros(K, numel(k)), 'hist', zeros(K, numel(k), max(opts.PhaseBins, 0)));
        for j = 1:numel(k)
            s = pairs(k(j), 1);  q = pairs(k(j), 2);
            [cs, cq, t] = i_onegrid(X, band, s, q, N, fs, cache);
            A  = abs(cq);
            ph = angle(cs);
            fi = min(floor(t / tile) + 1, K);
            e.vec(:, j) = accumarray(fi, A .* cos(ph), [K 1]) + 1i * accumarray(fi, A .* sin(ph), [K 1]);
            e.amp(:, j) = accumarray(fi, A, [K 1]);
            if opts.PhaseBins > 0
                bi = min(floor((ph + pi) / (2*pi) * opts.PhaseBins) + 1, opts.PhaseBins);
                e.hist(:, j, :) = reshape(accumarray([fi bi], A, [K opts.PhaseBins]), K, 1, []);
            end
        end
        P(end+1) = e;                                          %#ok<AGROW>
    end
end

% Both bands on ONE grid: the fast band's. The slow band's bins are placed in an array of
% the fast band's length, so the two come back at the same instants -- and each is
% re-modulated by its own first-bin frequency, undoing the sub-band demodulation.
function [cs, cq, t] = i_onegrid(X, band, s, q, N, fs, cache)
    nbq = band(q).n;
    Nm = tfb_smooth5(ceil(2 * nbq));
    t = (0:Nm-1)' * (N / Nm) / fs;
    cq = i_eval(X, band, q, Nm, N, t, cache);
    cs = i_eval(X, band, s, Nm, N, t, cache);
end

function c = i_eval(X, band, b, Nm, N, t, cache)
    key = sprintf('%d|%d', b, Nm);
    if isKey(cache, key), c = cache(key); return; end
    Y = zeros(Nm, 1);
    n = band(b).n;
    if n > Nm, error('ingest:couple:grid', 'Band %d does not fit the pair''s grid.', b); end
    Y(1:n) = X(band(b).i1:band(b).i2) .* band(b).H(:);
    c = ifft(Y) * (Nm / N);
    c = c .* exp(2i * pi * band(b).f0 * t);                    % undo the demodulation
    cache(key) = c;                                            %#ok<NASGU>
end

% A band's own band-pass: the sum of its members' squared gains, which is 1 inside the band
% and rolls off at the edges exactly as the tight frame does.
function band = i_bandbins(fb, bands, B, f)
    nB = height(bands);
    band = struct('i1', cell(1, nB), 'i2', [], 'n', [], 'H', [], 'f0', []);
    for b = 1:nB
        m = bands.scales{b};
        i1 = min(B(m, 1));  i2 = max(B(m, 2));
        H = zeros(i2 - i1 + 1, 1);
        for mm = m(:)'
            gm = gains(fb, mm, f(i1:i2));
            H = H + abs(gm(:)).^2;
        end
        band(b) = struct('i1', i1, 'i2', i2, 'n', i2 - i1 + 1, 'H', H, 'f0', f(i1));
    end
end

% Every ordered octave pair at least MinOctaves apart, slow first. ⚠ closer pairs are not
% carried by the fast band's rate and are refused rather than silently aliased.
function pairs = i_pairs(bands, spec, minOct)
    if isnumeric(spec)
        pairs = spec;
        if ~isempty(pairs)
            oc = log2(bands.fCenter(pairs(:,2)) ./ bands.fCenter(pairs(:,1)));
            bad = find(oc < minOct, 1);
            if ~isempty(bad)
                error('ingest:couple:pairs', ...
                      ['Bands %d and %d are %.2f octaves apart; a pair needs %g, because the ' ...
                       'fast band''s rate must carry twice the slow band''s top edge.'], ...
                      pairs(bad,1), pairs(bad,2), oc(bad), minOct);
            end
        end
        return
    end
    oct = find(strcmp(string(bands.kind), "octave"))';
    pairs = zeros(0, 2);
    for s = oct
        for q = oct
            if log2(bands.fCenter(q) / bands.fCenter(s)) >= minOct
                pairs(end+1, :) = [s q];                        %#ok<AGROW>
            end
        end
    end
end
% Author: Diellor Basha, 2026
