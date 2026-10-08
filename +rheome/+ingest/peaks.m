function P = peaks(x, bands, g, opts)
% INGEST.PEAKS  Per-tile spectral peak of each band, at the band's OWN level, by FFT.
%
%   P = rheome.ingest.peaks(x, bands, g)
%   P = rheome.ingest.peaks(x, bands, g, Pad=2, Window="hann")
%
% ⭐ WHY THIS IS NOT A ROLL-UP, AND DOES NOT NEED TO BE. A peak frequency is not mergeable:
% the loudest frequency of two joined tiles is not a function of each tile's peak. It does
% not matter, because the DIAGONAL gives every band exactly one home level -- the first tile
% at least its support long -- so a per-band measurement has one place to live and is
% computed there once. Not closed form, and not a problem.
%
% ⭐ WHY FFT RATHER THAN THE BANK. The frame bank already gives band ENERGY exactly, and its
% peak resolution is one member: at four voices that is a 1.19x step, far too coarse to call
% an alpha peak. A tile's own FFT resolves 1/T: a 2 s alpha tile gives 0.5 Hz bins, and a
% parabolic fit on the log-magnitude puts the peak between them. That is worth one more pass.
%
% ⚠ THE WINDOW COSTS AMPLITUDE, AND IT IS CORRECTED. A Hann window removes 50 % of the
% amplitude and leaks either side; the returned peak amplitude divides by the window's
% coherent gain, so it is the sinusoid's amplitude in the recording's units, not the
% spectrum's. Band POWER from this path is therefore NOT the store's `energy` -- that one is
% exact by Parseval over the whole record, this one is a windowed tile's own estimate, and
% the two agree to a few per cent on stationary content. Use `energy` for power; use this for
% where the power sits.
%
% INPUTS
%   x      [nT x 1] one channel, double
%   bands  the band table (rheome.ingest.bank)
%   g      the grid (rheome.ingest.grid)
%   Pad    zero-padding factor for the peak fit (2)
%   Window "hann" | "rect"
%
% OUTPUT (struct P), one entry per level that is some band's home:
%   .level      the level
%   .bands      [1 x nb] band ids whose natural level this is
%   .freq       [K x nb] peak frequency in Hz, NaN where the tile holds no samples
%   .amp        [K x nb] peak amplitude in the signal's units
%   .power      [K x nb] the band's power from this spectrum (units^2)
%
% See also: rheome.ingest.bank, rheome.ingest.build, rheome.select.derive
%
% Author: Diellor Basha, 2026

    arguments
        x double
        bands table
        g (1,1) struct
        opts.Pad    (1,1) double {mustBeGreaterThanOrEqual(opts.Pad, 1)} = 2
        opts.Window (1,1) string {mustBeMember(opts.Window, ["hann","rect"])} = "hann"
    end
    x = double(x(:));  nT = g.nT;  fs = g.fs;
    homes = unique(bands.naturalLevel(:))';
    P = struct('level', {}, 'bands', {}, 'freq', {}, 'amp', {}, 'power', {});

    for L = homes
        b = find(bands.naturalLevel == L)';
        K = g.K(L+1);  F = g.F * 2^L;
        nfft = max(8, 2^nextpow2(round(opts.Pad * F)));
        if opts.Window == "hann"
            w = 0.5 - 0.5*cos(2*pi*(0:F-1)'/F);                 % periodic Hann
        else
            w = ones(F, 1);
        end
        cg = sum(w) / F;                                        % coherent gain (amplitude)
        ng = sum(w.^2) / F;                                     % noise gain (power)
        fbin = (0:nfft-1)' * fs / nfft;

        % every tile of this level as a column, zero-padded at the record's end
        pad = K*F - nT;
        Xf = reshape([x; zeros(pad, 1)], F, K);
        nvalid = [repmat(F, 1, K - (pad > 0)), F - pad*(pad > 0)];
        Xf = Xf - mean(Xf, 1);                                  % a tile's DC is not its peak
        S = abs(fft(Xf .* w, nfft, 1)) / F;                     % [nfft x K]

        e = struct('level', L, 'bands', b, 'freq', nan(K, numel(b)), ...
                   'amp', nan(K, numel(b)), 'power', nan(K, numel(b)));
        for j = 1:numel(b)
            lo = bands.fLo(b(j));  hi = bands.fHi(b(j));
            k = find(fbin >= lo & fbin < hi & fbin <= fs/2);
            if isempty(k), continue; end
            Sb = S(k, :);
            [~, im] = max(Sb, [], 1);
            [fpk, apk] = i_interp(Sb, fbin(k), im);
            e.freq(:, j) = fpk(:);
            e.amp(:, j)  = 2 * apk(:) / cg;                     % one-sided, window-corrected
            % ⚠ PARSEVAL WITH ZERO PADDING. sum_k |X_k|^2 = nfft * sum_n |y_n|^2, so a
            % magnitude scaled by 1/F needs F/nfft to come back to mean power -- omitting it
            % overstates the band by exactly the padding factor (measured 3.45x at Pad=2,
            % against the bank's exact energy). The 2 is the one-sided fold, and 1/ng undoes
            % the window's power loss.
            e.power(:, j) = 2 * sum(Sb.^2, 1)' * F / (nfft * ng);
        end
        bad = nvalid(:) < F/2;                                  % a tile the record barely reaches
        e.freq(bad, :) = NaN;  e.amp(bad, :) = NaN;  e.power(bad, :) = NaN;
        P(end+1) = e;                                           %#ok<AGROW>
    end
end

% ⚠ A PEAK BETWEEN BINS. The maximum bin is not the peak: a parabola through the log
% magnitudes of the bin and its neighbours puts it where it belongs, to a fraction of a bin,
% which is what makes a 0.5 Hz grid usable for an alpha peak. At an edge bin there is no
% neighbour, so the bin itself is the answer.
function [f, a] = i_interp(S, fb, im)
    K = size(S, 2);  n = size(S, 1);
    f = zeros(1, K);  a = zeros(1, K);
    for c = 1:K
        i = im(c);
        a(c) = S(i, c);  f(c) = fb(i);
        if i > 1 && i < n
            y = log(max(S(i-1:i+1, c), realmin));
            den = y(1) - 2*y(2) + y(3);
            if den ~= 0
                d = 0.5 * (y(1) - y(3)) / den;                  % in bins, |d| <= 0.5
                d = max(min(d, 0.5), -0.5);
                f(c) = fb(i) + d * (fb(2) - fb(1));
                a(c) = exp(y(2) - 0.25 * (y(1) - y(3)) * d);
            end
        end
    end
end
% Author: Diellor Basha, 2026
