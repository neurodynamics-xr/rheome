function T0 = reduce(X, fb, bands, g, frame, opts)
% INGEST.REDUCE  Level-0 tile statistics from samples and from CWT coefficients, under a mask.
%
%   T0 = rheome.ingest.reduce(X, fb, bands, g, frame)
%   T0 = rheome.ingest.reduce(X, fb, bands, g, frame, Mask=m, Precision="single", MaxBytes=2e9)
%
% The only function that touches samples. Per channel: the whole record goes through the
% bank ONCE (wt), per-sample band series are formed from the band's scales, and every
% series is reshaped into level-0 frames and reduced. With a timefilterbank (the default)
% each member arrives at its OWN rate
% and is framed by its instants: energy is sum |coef|^2 * weight per frame (Parseval),
% envMax the max over the member's instants; a masked instant is one whose nearest
% record sample is masked. 7x faster than the full-rate path, no paging. Every output is a sum, a count, a
% max or a min over samples, so rheome.ingest.rollup can form every coarser level from this one
% without the signal (design §2, §3.1).
%
% ⚠ THE MASK ACTS HERE, NOT BEFORE THE TRANSFORM. The CWT sees the record as stored;
% masked samples are then excluded from every sum and count. Zeroing samples before the
% transform would inject a step at every mask edge into every band (design §3.7).
% Coefficients near a masked run are therefore computed from the unmasked signal.
%
% ⚠ Statistics that depend only on the mask and the grid are stored ONCE per record, not
% per channel: .n and .nCoi are [K0 x 1] and [K0 x nB]. They gain a channel axis the day a
% per-channel mask exists.
%
% ⚠ Hand-rolled on purpose. Measured on 600 s x 270 channels at level 0:
% signalTimeFeatureExtractor 17.1 s, this reshape-and-reduce 0.08 s, identical results;
% and the extractor rejects any non-finite sample, so it cannot honour a mask. It remains
% the ORACLE in tests/tIngestReduce.m (inventory §3, design §3.3).
%
% Cone of influence: analytic, from the record's edges (rheome.ingest.cone). Band j is inside the
% cone at t iff its lowest scale is below coi(t) = min(fmax, coiConst/dist(t)). nCoi counts
% the valid samples of a frame for which that holds. wt's own coi agrees with the analytic
% one to the fit's 3e-4; the analytic form is used so that the paged path counts the same.
%
% INPUTS:
%   X      [nT x C] samples (any float class; the transform runs at Precision)
%   fb     cwtfilterbank from rheome.ingest.bank (SignalLength must equal nT)
%   bands  band table from rheome.ingest.bank
%   g      grid from rheome.ingest.grid
%   frame  frame constants from rheome.ingest.bank (for the cone)
%   Mask   [nT x 1] logical, true = valid sample (default all true)
%   Precision  "single" | "double" for the CWT buffer (sums are always double)
%   MaxBytes   byte guard on the per-channel CWT buffer -> error ingest:reduce:bytes
% OUTPUT (struct T0, level 0, K0 frames):
%   .n      [K0 x 1] uint32         valid samples per frame
%   .sumX   [K0 x C] double         sum x
%   .sumX2  [K0 x C] double         sum x^2
%   .absMax [K0 x C] single         max |x|        (0 on an empty frame)
%   .min    [K0 x C] single         min x          (+Inf on an empty frame)
%   .max    [K0 x C] single         max x          (-Inf on an empty frame)
%   .energy [K0 x C x nB] double    sum_t sum_{m in band} |W_m(t)|^2
%   .envMax [K0 x C x nB] single    max_t max_{m in band} |W_m(t)|
%   .nCoi   [K0 x nB] uint32        valid samples inside the cone of influence
%   .level = 0, .K = K0
%
% See also: rheome.ingest.bank, rheome.ingest.grid, rheome.ingest.rollup, cwtfilterbank/wt
%
% Author: Diellor Basha, 2026

    arguments
        X      {mustBeNumeric}
        fb     (1,1)                       % cwtfilterbank (full rate) or timefilterbank (sub-band)
        bands  table
        g      (1,1) struct
        frame  (1,1) struct
        opts.Mask      logical = logical([])
        opts.Precision (1,1) string {mustBeMember(opts.Precision, ["single","double"])} = "single"
        opts.MaxBytes  (1,1) double {mustBePositive} = 2e9
    end

    [nT, C] = size(X);
    if nT ~= g.nT
        error('ingest:reduce:length', 'X has %d samples but the grid was built for %d.', nT, g.nT);
    end
    mask = opts.Mask;
    if isempty(mask), mask = true(nT, 1); end
    mask = mask(:);
    if numel(mask) ~= nT
        error('ingest:reduce:mask', 'Mask has %d elements for %d samples.', numel(mask), nT);
    end

    fc = centerFrequencies(fb);  fc = fc(:)';
    nF = numel(fc);
    isFrame = isa(fb, 'rheome.timefilterbank');
    bytesPer = 8 + 8 * (opts.Precision == "double");        % complex
    bytes = 5 * nF * nT * bytesPer;                          % wt's working set, measured ~5x the buffer
    if ~isFrame && bytes > opts.MaxBytes
        error('ingest:reduce:bytes', ...
              'The per-channel CWT working set needs %.3g bytes (5 x %d filters x %d samples, %s); MaxBytes is %.3g.', ...
              bytes, nF, nT, opts.Precision, opts.MaxBytes);
    end

    F  = g.F;  K0 = g.K0;  nB = height(bands);
    pad = K0*F - nT;
    frames = @(v, fill) reshape([double(v(:)); repmat(fill, pad, 1)], F, K0);   % [F x K0]

    mk = frames(mask, 0) > 0;                                  % valid-sample mask per frame
    n  = uint32(sum(mk, 1))';

    T0 = struct();
    T0.n      = n;
    T0.sumX   = zeros(K0, C);
    T0.sumX2  = zeros(K0, C);
    T0.sumAbs  = zeros(K0, C);
    T0.sumSqrt = zeros(K0, C);
    T0.sumX3   = zeros(K0, C);
    T0.sumX4   = zeros(K0, C);
    T0.absMax = zeros(K0, C, 'single');
    T0.min    = zeros(K0, C, 'single');
    T0.max    = zeros(K0, C, 'single');
    T0.energy = zeros(K0, C, nB);
    T0.envMax = zeros(K0, C, nB, 'single');
    T0.nCoi   = zeros(K0, nB, 'uint32');
    T0.level  = 0;
    T0.K      = K0;

    lowest = cellfun(@(s) fc(s(end)), bands.scales);            % lowest scale centre per band
    valid  = mask';                                             % [1 x nT]
    inside = rheome.ingest.cone(nT, g.fs, frame, lowest);                % [nB x nT]
    for b = 1:nB
        T0.nCoi(:, b) = uint32(sum(frames(inside(b, :) & valid, 0), 1))';
    end

    for c = 1:C
        x = double(X(:, c));

        % all-pass, from the samples
        xm = x;  xm(~mask) = 0;
        Xf = frames(xm, 0);
        T0.sumX(:, c)   = sum(Xf, 1)';
        T0.sumX2(:, c)  = sum(Xf.^2, 1)';
        % ⭐ FOUR MORE SUMS, AND THEY MERGE LIKE THE OTHERS. Between them they carry the whole
        % shape family a time-domain extractor reports -- skewness and kurtosis from the third
        % and fourth moments, shape/impulse/clearance factors from sum|x| and sum sqrt|x| --
        % none of which can be recovered from mean and rms alone, and all of which roll up the
        % pyramid exactly because a sum of sums is a sum.
        T0.sumAbs(:, c)  = sum(abs(Xf), 1)';
        T0.sumSqrt(:, c) = sum(sqrt(abs(Xf)), 1)';
        T0.sumX3(:, c)   = sum(Xf.^3, 1)';
        T0.sumX4(:, c)   = sum(Xf.^4, 1)';
        T0.absMax(:, c) = ing_bound(max(abs(Xf), [], 1), 'up')';        % a bound, not a round
        xm = x;  xm(~mask) = +Inf;   T0.min(:, c) = ing_bound(min(frames(xm, +Inf), [], 1), 'down')';
        xm = x;  xm(~mask) = -Inf;   T0.max(:, c) = ing_bound(max(frames(xm, -Inf), [], 1), 'up')';

        % per band, from the coefficients of the record as stored
        if isFrame
            % sub-band path: each member at its own rate; frames by the member's instants
            Cc = wt(fb, x);
            for b = 1:nB
                for mIdx = bands.scales{b}
                    cm = Cc(mIdx);
                    mag = abs(cm.coef);
                    si  = min(round(cm.t * g.fs) + 1, nT);          % nearest record sample
                    ok  = mask(si);
                    fi  = min(floor(cm.t / g.FrameFloor) + 1, K0);
                    e = mag.^2 * cm.weight;  e(~ok) = 0;
                    v = mag;                 v(~ok) = 0;
                    T0.energy(:, c, b) = T0.energy(:, c, b) + accumarray(fi, e, [K0 1]);
                    T0.envMax(:, c, b) = max(T0.envMax(:, c, b), ing_bound(accumarray(fi, v, [K0 1], @max), 'up'));
                end
            end
            continue
        end
        if opts.Precision == "single", xw = single(x); else, xw = x; end
        W = wt(fb, xw);                                         % [nF x nT]
        for b = 1:nB
            Wb = abs(W(bands.scales{b}, :));                    % [nS x nT]
            e  = sum(double(Wb).^2, 1);   e(~valid) = 0;
            v  = max(Wb, [], 1);          v(~valid) = 0;
            T0.energy(:, c, b) = sum(frames(e, 0), 1)';
            T0.envMax(:, c, b) = ing_bound(max(frames(v, 0), [], 1), 'up')';
        end
    end
end

% Author: Diellor Basha, 2026
