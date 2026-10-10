function T = wavelettile(db, varargin)
% SELECTION.WAVELETTILE  Wavelets as a tiling: centre and time SUPPORT where a tile has centre and span.
%
%   T = rheome.selection.wavelettile(db)                       % every band the store carries
%   T = rheome.selection.wavelettile(db, Bands=6)              % alpha only
%   T = rheome.selection.wavelettile(db, Energy=0.99, Positions=true)
%
% ⭐ THE IDEA. A dyadic tile is a hard window: a centre and a span, and a set of them partitions the
% line. A wavelet is the same pair of numbers -- a centre and a time SUPPORT -- with the span replaced
% by a decay. So a band can be selected by the wavelet that best supports it, and that wavelet's
% support localises it in time: one object navigates frequency and time together, which is what a tile
% and a band do only when paired by hand.
%
% ⭐⭐ THE LEVEL STRUCTURE IS IDENTICAL TO THE HARD TILING, AND THAT IS THE WHOLE POINT. One wavelet's
% support is EXACTLY two of the next level's: supportCycles is scale-invariant (15.051) and fc halves
% per octave, so supportSec goes 0.166 / 0.333 / 0.665 / 1.330 / 2.661 / 5.322 s, each exactly twice the
% last. A wavelet at level L sits above two at L+1, four at L+2, down to the sampling rate -- the same
% nesting a dyadic tile has. `nested = true` in rheome.selection.registry means that, and nothing looser.
%
% ⭐ THERE ARE TWO STRIDES, FOR TWO PURPOSES, and conflating them is what made the first version of
% this function confusing:
%   strideSelect = supportSec        redundancy 1     the TILE LATTICE with a smooth window instead of
%                                                     a hard one. This is the navigation stride: one
%                                                     position per support, levels nesting two-to-one.
%   strideFrame  = 1/rate            redundancy 12.5  the bank's own translation step, set by the
%                                                     BANDWIDTH of the coefficient sequence. Required
%                                                     for the frame to invert.
% ⚠ The selection lattice is a sub-lattice of the frame lattice -- roughly every 12th position -- so it
% indexes and navigates but does NOT reconstruct: dropping to it loses the frame property. Select on
% strideSelect; reconstruct through the full bank.
%
% ⚠ AND THE TWO STRIDES DO NOT AGREE ON WHAT "ONE PER TILE" MEANS. The dyadic tile is 1.50 supports
% long (supportPerTile = 0.665), because the tile is the support rounded UP to the next power of two.
% So a wavelet at each tile centre covers two thirds of its tile and a covering needs 1.5 per tile;
% `positionsPerTileSelect` says so rather than leaving it to be discovered.
%
% ⚠ "stride from the bandwidth" and "stride from the dyadic level" are the SAME statement here, not two
% accounts: constant Q makes bandwidth proportional to fc, so a bandwidth-derived stride is
% automatically dyadic. Neither derivation is more fundamental.
%
% ⚠⚠ A WAVELET TILING DOES NOT PARTITION, so nothing summed over it is exact. `freq_band` and
% `time_tile` merge -- a parent is exactly its children -- and this does not: measured, one voice holds
% 0.245 of its band's response and a band's own four voices hold only 0.835 of the bank's response over
% that band, the rest being neighbouring voices leaking in. Use this to SELECT and to navigate; use the
% tiles when a total has to add up. rheome.selection.registry marks it kind='kernel', merges=false.
%
% ⚠ THE SUPPORT IS AN ENERGY THRESHOLD, NOT AN ANALYTIC WIDTH, and the two do not agree numerically.
% `support(tfb)` is the fraction-of-energy definition; `widths` on the graph side is analytic
% (sigma = sqrt(2t), exact). Report which one a number came from -- this table's `Energy` column says.
% ⚠ I ASSERTED THE SUPPORT DEPENDS ON SignalLength AND IT DOES NOT -- measured 15.0514 cycles at
% 2^14, 2^16 and 2^18, identical to every digit. That warning came from `cwtfilterbank/waveletsupport`,
% which is a different function with a different contract; `timefilterbank/support` returns the
% scale-invariant constant and is length-independent. The column is kept as provenance, not as a hazard.
%
% ⚠⚠ WHAT THE SUPPORT REALLY DEPENDS ON IS THE ENERGY THRESHOLD, and strongly: measured 5.92 / 8.31 /
% 15.05 cycles at Energy = 0.95 / 0.99 / 0.999, so `supportPerTile` swings 0.26 -> 0.37 -> 0.67. Any
% statement of the form "the wavelet supporting alpha is 1.33 s long" is a statement about 0.999 of its
% energy and nothing else. Quote the threshold with the number.
%
% ⭐ THE EXCHANGE RATE WITH THE HARD TILING, measured: a wavelet spans 15.05 cycles of its own centre
% frequency (scale-invariant) and a constant-Q tile spans 22.63, so A TILE IS 1.50 WAVELET SUPPORTS
% LONG -- the tile is the support rounded up to the next dyadic level. `supportPerTile` is that ratio
% per row and should read 0.67 in every octave.
%
% NAME-VALUE
%   Bands      store band ids to include (default: all with a natural level in the grid)
%   Energy     fraction of |psi|^2 the support must hold (default 0.999)
%   Positions  also return the lattice of centre times per band (default false)
%   Stride     "select" (default, = support, the tile lattice) | "frame" (= 1/rate, the bank's own)
%   SignalLength  samples the support is measured over (default 2^16) -- see the note above
%
% OUTPUT (one row per band, or per (band, position) with Positions=true)
%   band, fLo, fHi, fCenter, q, supportSec, supportCycles, strideSelect, strideFrame,
%   redundancyFrame, positionsPerTileSelect, positionsPerTileFrame, childrenPerSupport,
%   dyadicLevel, tileSec, cyclesPerTile, supportPerTile, energy, signalLength, tCenter (Positions)
%
% See also: rheome.selection.registry, rheome.select.ladder, rheome.geom.ladder, rheome.timefilterbank/support,
%           @graphfilterbank/widths
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Bands', [], @isnumeric);
    p.addParameter('Energy', 0.999, @(x) isscalar(x) && x > 0 && x < 1);
    p.addParameter('Positions', false, @islogical);
    p.addParameter('Stride', "select", @(x) any(strcmpi(x, ["select" "frame"])));
    p.addParameter('SignalLength', 2^16, @isscalar);
    p.parse(varargin{:});
    o = p.Results;

    L   = rheome.select.ladder(db);
    fs  = db.meta.fs;
    tfb = rheome.timefilterbank(o.SignalLength, SamplingFrequency=fs, VoicesPerOctave=4, Anchor=1);
    fc  = centerFrequencies(tfb);
    sup = support(tfb, 'Energy', o.Energy);               % ⚠ a STRUCT; .timeTimesFc is the invariant
    cyc = sup.timeTimesFc;

    rows = {};
    for i = 1:height(L)
        b = L.band_id(i);
        if ~isempty(o.Bands) && ~ismember(b, o.Bands), continue; end
        if L.natural_level(i) < 0 || ~isfinite(L.tile_s(i)), continue; end
        % the voice that best supports this band: the one nearest its geometric centre
        v = find(fc >= L.f_lo(i) & fc < L.f_hi(i));
        if isempty(v), continue; end
        [~, ic] = min(abs(log(fc(v)/L.f_center(i))));  mc = v(ic);
        supS = cyc / fc(mc);                               % support in seconds, scale-invariant form
        sSel = supS;                                       % ⭐ navigation: one position per support
        sFrm = 1/L.rate_hz(i);                             % the bank's own step, set by bandwidth
        rows{end+1} = table(b, L.f_lo(i), L.f_hi(i), fc(mc), L.q(i), supS, cyc, sSel, sFrm, ...
            supS/sFrm, L.tile_s(i)/sSel, L.tile_s(i)/sFrm, ...
            L.natural_level(i), L.tile_s(i), L.cycles_per_tile(i), supS/L.tile_s(i), ...
            o.Energy, o.SignalLength, ...
            'VariableNames', {'band','fLo','fHi','fCenter','q','supportSec','supportCycles', ...
            'strideSelect','strideFrame','redundancyFrame','positionsPerTileSelect', ...
            'positionsPerTileFrame','dyadicLevel','tileSec','cyclesPerTile','supportPerTile', ...
            'energy','signalLength'}); %#ok<AGROW>
    end
    if isempty(rows), error('selection:wavelettile:none', 'No band matched.'); end
    T = vertcat(rows{:});
    % ⭐ the nesting, measured rather than asserted: how many of the NEXT LEVEL's supports fit in this
    % one. Must be exactly 2 between consecutive octave bands.
    % ⚠ BAND IDS RISE AS FREQUENCY FALLS on this bank (band 1 is 215-300 Hz, band 11 is 0.25-0.5 Hz),
    % so the child -- finer in time, shorter support -- is the band one LOWER in id. Computing
    % support(i)/support(i+1) gives 0.5 and looks like a failed check; it is the same fact upside down.
    T.childrenPerSupport = nan(height(T),1);
    for i = 2:height(T)
        if T.band(i) == T.band(i-1) + 1
            T.childrenPerSupport(i) = T.supportSec(i) / T.supportSec(i-1);
        end
    end

    if o.Positions
        % ⭐ the lattice: centre times at the band's own stride, which is what makes it a TILING
        out = {};
        nT = db.grid.nT / fs;
        useSel = strcmpi(o.Stride, "select");
        for i = 1:height(T)
            st = T.strideSelect(i);  if ~useSel, st = T.strideFrame(i); end
            tc = (st/2 : st : nT)';
            r  = repmat(T(i,:), numel(tc), 1);
            r.tCenter = tc;
            out{end+1} = r; %#ok<AGROW>
        end
        T = vertcat(out{:});
    end
end

% Author: Diellor Basha, 2026
