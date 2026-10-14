function P = eventlocked(X, fs, anchors, filt, opts)
% DETECT.EVENTLOCKED  Time-frequency power around event anchors (onsets, offsets, random times).
%
%   P = rheome.detect.eventlocked(X, fs, anchors, filt)
%   P = rheome.detect.eventlocked(X, fs, anchors, filt, Window=[-1 1], Freqs=1:45, Cycles=5, Rate=100)
%       X        [C x nS] data (sensors)
%       anchors  [nE x 1] times (s, sample 1 = 0 s), e.g. the events' tOn
%       filt     {nE x 1} cell of [k x C] spatial filters, or one [k x C] for all (e.g. eye(C): every channel);
%                the power is summed over the filter's k outputs -- for the event-tile source, the minimum-norm
%                kernel rows at the event's tile centres, each tile's 3 rows scaled by sqrt(its area share)
%
% Complex Morlet wavelets (Cycles cycles, Gaussian in frequency, sigma_f = f / Cycles) on a segment padded
% by 3 sigma_t of the lowest frequency on either side, so the window carries no edge.
% P [nF x nT x nE] power; an anchor whose padded segment leaves the data is all NaN. Times P.t below.
% Use: average P(:,:,e) ./ B(:,e) over events, B the same event's power at matched random anchors averaged
% over time and anchors (rheome.scale.measure_events), in dB.
%
% OUTPUT struct P: .P [nF x nT x nE]   .f [nF x 1]   .t [1 x nT] s from the anchor
%
% See also: rheome.detect.eventtable, rheome.scale.measure_events, rheome.scale.measure_slowosc
%
% Author: Diellor Basha, 2026

    arguments
        X double
        fs (1,1) double {mustBePositive}
        anchors (:,1) double
        filt
        opts.Window (1,2) double = [-1 1]
        opts.Freqs (:,1) double {mustBePositive} = (1:45)'
        opts.Cycles (1,1) double {mustBePositive} = 5
        opts.Rate (1,1) double {mustBePositive} = 100
    end
    f = opts.Freqs;  nE = numel(anchors);  nS = size(X, 2);
    pad = 3 * opts.Cycles / (2*pi*min(f));
    i0 = round((opts.Window(1) - pad) * fs);  n = round((diff(opts.Window) + 2*pad) * fs) + 1;
    fk = (0:n-1)' * fs / n;
    Psi = 2 * exp(-0.5 * ((fk - f') ./ (f' / opts.Cycles)).^2);         % [n x nF], analytic
    keep = round((pad + (0:1/opts.Rate:diff(opts.Window))) * fs) + 1;   % window samples, at Rate
    P.t = opts.Window(1) + (keep - keep(1)) / fs;  P.f = f;
    P.P = nan(numel(f), numel(keep), nE);
    for e = 1:nE
        s = round(anchors(e) * fs) + 1 + i0 + (0:n-1);
        if s(1) < 1 || s(end) > nS || ~isfinite(anchors(e)), continue, end
        W = filt;  if iscell(filt), W = filt{e}; end
        x = (W * X(:, s)).';                                            % [n x k]
        Z = fft(x - mean(x, 1), [], 1);
        acc = zeros(numel(keep), numel(f));
        for j = 1:size(Z, 2)
            c = ifft(Z(:, j) .* Psi, [], 1);  acc = acc + abs(c(keep, :)).^2;
        end
        P.P(:, :, e) = acc.';
    end
end

% Author: Diellor Basha, 2026
