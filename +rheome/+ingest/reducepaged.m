function T0 = reducepaged(readfcn, C, bands, g, frame, plan, opts)
% INGEST.REDUCEPAGED  Level-0 statistics computed page by page, one sub-bank per job.
%
%   T0 = rheome.ingest.reducepaged(readfcn, C, bands, g, frame, plan)
%   T0 = rheome.ingest.reducepaged(..., Mask=m, Precision="single", MaxBytes=2e9)
%
% The paged twin of rheome.ingest.reduce: same output, same definitions, but the record is never
% transformed whole. For each job of the plan (rheome.ingest.pages) and each page at the job's
% level, a span of core +/- halo is read through readfcn, a SUB-BANK anchored on the
% job's top master scale is applied (its scales are the master's to 4e-16), and the
% band series are reduced over the CORE only -- the halo belongs to the neighbours and is
% never counted (the same core/span rule @pagedrecording pins). Rows of the level-0
% arrays are written in place, so the result concatenates across pages without merging.
%
% ⚠ Paged coefficients differ from whole-record ones by the wavelet tail outside the
% halo. waveletsupport's support is a threshold, not a zero. Measured against the
% whole-record store. The page policy is therefore
% in the config hash.
%
% ⚠ nCoi is analytic and identical to rheome.ingest.reduce's (rheome.ingest.cone): a page's own edges
% are covered by its halo and are NOT a cone. Where the halo is clipped by the record
% edge, the cone count is what says so.
%
% INPUTS:
%   readfcn  @(a, b) -> [b-a+1 x C] samples a..b (1-based, inclusive), any float class
%   C        channels readfcn returns
%   bands, g, frame   from rheome.ingest.bank / rheome.ingest.grid
%   plan     from rheome.ingest.pages
%   Mask, Precision, MaxBytes   as in rheome.ingest.reduce
% OUTPUT:
%   T0  as rheome.ingest.reduce
%
% See also: rheome.ingest.reduce, rheome.ingest.pages, rheome.ingest.bank, rheome.ingest.cone
%
% Author: Diellor Basha, 2026

    arguments
        readfcn (1,1) function_handle
        C       (1,1) double {mustBeInteger, mustBePositive}
        bands   table
        g       (1,1) struct
        frame   (1,1) struct
        plan    table
        opts.Mask      logical = logical([])
        opts.Precision (1,1) string {mustBeMember(opts.Precision, ["single","double"])} = "single"
        opts.MaxBytes  (1,1) double {mustBePositive} = 2e9
    end

    if isfield(frame, 'bank') && strcmp(frame.bank, 'frame')
        error('ingest:reducepaged:bank', 'The paged path is for the full-rate Morse bank; a rheome.timefilterbank needs no paging (use rheome.ingest.reduce).');
    end
    nT = g.nT;  fs = g.fs;  F = g.F;  K0 = g.K0;  nB = height(bands);
    mask = opts.Mask;
    if isempty(mask), mask = true(nT, 1); end
    mask = mask(:);
    if numel(mask) ~= nT
        error('ingest:reducepaged:mask', 'Mask has %d elements for %d samples.', numel(mask), nT);
    end
    if any(plan.bufferBytes > opts.MaxBytes)
        error('ingest:reducepaged:bytes', 'A job needs %.3g bytes per channel; MaxBytes is %.3g.', ...
              max(plan.bufferBytes), opts.MaxBytes);
    end
    fc = frame.fc;  V = numel(bands.scales{1});
    lowest = cellfun(@(s) fc(s(end)), bands.scales);

    T0 = struct();
    T0.n      = uint32(i_framesum(double(mask), F, K0, 0));
    T0.sumX   = zeros(K0, C);
    T0.sumX2  = zeros(K0, C);
    T0.absMax = zeros(K0, C, 'single');
    T0.min    = zeros(K0, C, 'single');
    T0.max    = zeros(K0, C, 'single');
    T0.energy = zeros(K0, C, nB);
    T0.envMax = zeros(K0, C, nB, 'single');
    T0.nCoi   = zeros(K0, nB, 'uint32');
    T0.level  = 0;
    T0.K      = K0;

    inside = rheome.ingest.cone(nT, fs, frame, lowest);
    for b = 1:nB
        T0.nCoi(:, b) = uint32(i_framesum(double(inside(b, :)' & mask), F, K0, 0));
    end

    % all-pass statistics: read the record once in cores at the lowest job level
    Lap = min(plan.level);
    for p = 1:g.K(Lap+1)
        [a, b] = i_core(p, Lap, g);
        X = double(readfcn(a, b));  m = mask(a:b);
        k = i_frames_of(a, b, F);
        Xm = X;  Xm(~m, :) = 0;
        T0.sumX(k, :)   = i_red(Xm, F, @sum, 0);
        T0.sumX2(k, :)  = i_red(Xm.^2, F, @sum, 0);
        % ⚠ SAME OUTWARD ROUNDING AS THE WHOLE-RECORD PATH (rheome.ingest.bound). Rounding to
        % nearest here and outward there would put the two paths one ulp apart, which is
        % exactly what the paged-vs-whole identity test measures.
        T0.absMax(k, :) = ing_bound(i_red(abs(Xm), F, @max, 0), 'up');
        Xm = X;  Xm(~m, :) = +Inf;  T0.min(k, :) = ing_bound(i_red(Xm, F, @min, +Inf), 'down');
        Xm = X;  Xm(~m, :) = -Inf;  T0.max(k, :) = ing_bound(i_red(Xm, F, @max, -Inf), 'up');
    end

    % band statistics: per job, per page, per channel
    half = 2^(1/(2*V));
    for r = 1:height(plan)
        L = plan.level(r);  halo = plan.haloSamples(r);
        sc = plan.scales(r, 1):plan.scales(r, 2);          % master indices, top first
        limits = [fc(sc(end)) / half * 0.999, fc(sc(1))];
        cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
        jb = plan.bands{r};
        rel = cellfun(@(s) s - sc(1) + 1, bands.scales(jb), 'UniformOutput', false);   % rows in W
        for p = 1:plan.nPages(r)
            [a, b] = i_core(p, L, g);
            s1 = max(1, a - halo);  s2 = min(nT, b + halo);
            nSpan = s2 - s1 + 1;
            if ~isKey(cache, nSpan)
                cache(nSpan) = i_subbank(nSpan, fs, limits, fc(sc), frame, V);
            end
            sub = cache(nSpan);
            X = readfcn(s1, s2);
            if opts.Precision == "single", X = single(X); else, X = double(X); end
            ci = (a - s1 + 1):(b - s1 + 1);                  % core columns within the span
            m  = mask(a:b)';
            k  = i_frames_of(a, b, F);
            for c = 1:C
                W = wt(sub, X(:, c));                        % [nS x nSpan]
                for i = 1:numel(jb)
                    Wb = abs(W(rel{i}, ci));
                    e = sum(double(Wb).^2, 1);  e(~m) = 0;
                    v = max(Wb, [], 1);         v(~m) = 0;
                    T0.energy(k, c, jb(i)) = i_red(e', F, @sum, 0);
                    T0.envMax(k, c, jb(i)) = ing_bound(i_red(v', F, @max, 0), 'up');
                end
            end
        end
    end
end

function [a, b] = i_core(p, L, g)
    n = g.F * 2^L;
    a = (p-1) * n + 1;  b = min(p * n, g.nT);
end

function k = i_frames_of(a, b, F)
    k = ((a-1)/F + 1):ceil(b / F);
end

function R = i_red(X, F, op, fill)
% reduce [n x C] into frames of F rows (last frame padded with fill) -> [nfr x C]
    [n, C] = size(X);
    nfr = ceil(n / F);
    if nfr * F > n, X(n+1:nfr*F, :) = fill; end
    Y = reshape(X, F, nfr, C);
    if isequal(op, @sum), R = sum(Y, 1); else, R = op(Y, [], 1); end
    R = reshape(R, nfr, C);
end

function s = i_framesum(v, F, K0, fill)
    n = numel(v);  v = [v(:); repmat(fill, K0*F - n, 1)];
    s = sum(reshape(v, F, K0), 1)';
end

function sub = i_subbank(nSpan, fs, limits, fcWanted, frame, V)
    sub = cwtfilterbank('SignalLength', nSpan, 'SamplingFrequency', fs, ...
                        'Wavelet', frame.wavelet, 'VoicesPerOctave', V, 'FrequencyLimits', limits);
    got = centerFrequencies(sub);  got = got(:)';
    nW = numel(fcWanted);
    if numel(got) < nW || max(abs(got(1:nW) - fcWanted) ./ fcWanted) > 1e-9
        error('ingest:reducepaged:subbank', ...
              'The sub-bank for span %d does not reproduce master scales %.4g..%.4g Hz.', nSpan, fcWanted(1), fcWanted(end));
    end
    if numel(got) > nW
        % an extra scale below the wanted ones is harmless: it is never read
    end
end
% Author: Diellor Basha, 2026
