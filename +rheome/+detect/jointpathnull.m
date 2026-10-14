function N = jointpathnull(Ysur, G, nF, rate, opts)
% DETECT.JOINTPATHNULL  Normalisation, baseline, threshold, cross-band criterion and held-out false rates for rheome.detect.jointpath.
%
%   N = rheome.detect.jointpathnull(Ysur, G, nF, rate, LogF=log2fc)              % Alpha 0.1/s, BaseQ 0.9, CrossQ 0.95
%   N = rheome.detect.jointpathnull(Ysur, G, nF, rate, LogF=log2fc, Alpha=0.05, CrossOct=1)   % a fixed criterion
%
% Ysur is a NULL band-envelope field -- a surrogate of the data that destroys the alignment of places and KEEPS
% the coupling of bands at one place (one circular shift per place, shared by every band: overlapping bands
% share envelope variance, and a per-band shift would remove it and under-state the cross-band false rate) --
% cut into windows of nF frames and concatenated: [nTile x nBand x nF*nWin]. As rheome.detect.tilepathnull:
% odd windows CALIBRATE, even windows TEST.
%   norm       [nTile x nBand] median of each series over the calibration windows (jointpath's Norm)
%   baseline   the BaseQ quantile of the normalised calibration values
%   threshold  the best-path score at most Alpha*WindowS calibration windows exceed (Alpha false paths/s)
%   crossOct   the cross-band criterion on |netOct| (jointpath): the CrossQ quantile of |netOct| over the
%              PASSING calibration paths, unless CrossOct is given
%   falseTestS passing best paths per second on the TEST windows (any path); falseTestW, per window
%   crossTestS the same, counting only paths with |netOct| >= crossOct: the false rate of the claim this method
%              exists to make, "a pattern moved across frequency"; crossTestW, per window
%   bandOcc    [1 x nBand] fraction of best-path frames in each band, over all null windows. ⚠ Low bands have
%              slower envelopes on a common frame grid, so noise paths sit mostly there: the threshold is set
%              by low-band paths and the sensitivity differs by band.
%
% OUTPUT (struct N): .norm .baseline .threshold .crossOct .crossQ .falseTestS .falseTestW .crossTestS .crossTestW
%   .bandOcc .alpha .windowS .calScores .calNet (|netOct| of passing calibration paths) .testScores .testSpan
%   .testSweep .testNet [nWin/2 x 1]
%
% See also: rheome.detect.jointpath, rheome.detect.tilepathnull
%
% Author: Diellor Basha, 2026

    arguments
        Ysur double
        G struct
        nF (1,1) double {mustBeInteger, mustBePositive}
        rate (1,1) double {mustBePositive}
        opts.LogF double = []
        opts.Alpha (1,1) double {mustBePositive} = 0.1
        opts.BaseQ (1,1) double {mustBeInRange(opts.BaseQ, 0, 1)} = 0.9
        opts.WindowS = []
        opts.CrossOct double = []
        opts.CrossQ (1,1) double {mustBeInRange(opts.CrossQ, 0, 1)} = 0.95
    end
    [nTile, nB, nT] = size(Ysur);  nWin = floor(nT/nF);
    if nWin < 4, error('detect:jointpathnull:short', 'Need >= 4 null windows, got %d.', nWin); end
    W = opts.WindowS;  if isempty(W), W = nF/rate; end
    cal = 1:2:nWin;  tst = 2:2:nWin;
    ci = reshape((cal(:) - 1)*nF + (1:nF), 1, []);
    nrm = median(Ysur(:, :, ci), 3);  nrm(~(nrm > 0)) = realmin;      % a dead series must not divide by 0
    Zc = Ysur(:, :, ci) ./ nrm;  b = quantile(Zc(:), opts.BaseQ);
    sc = zeros(nWin, 1);  sp = nan(nWin, 1);  sw = sp;  nt = sp;  occ = zeros(1, nB);
    for w = 1:nWin
        p = rheome.detect.jointpath(Ysur(:, :, (w-1)*nF + (1:nF)), G, Norm=nrm, Baseline=b, LogF=opts.LogF, SampleRate=rate);
        if p.nTracks
            x = p.tracks(1);  sc(w) = x.score;  sp(w) = x.octSpan;  sw(w) = x.sweepOctS;  nt(w) = abs(x.netOct);
            occ = occ + accumarray(x.bands(:), 1, [nB 1])';
        end
    end
    s = sort(sc(cal), 'descend');  k = floor(opts.Alpha*W*numel(s));
    if k >= numel(s), thr = 0; else, thr = s(k + 1); end
    pc = sc(cal) > thr;  cn = nt(cal(pc));
    xo = opts.CrossOct;
    if isempty(xo)                                                   % calibrated: no passing path -> no cross claim
        if isempty(cn), xo = Inf; else, xo = quantile(cn, opts.CrossQ); end
    end
    ps = sc(tst) > thr;  px = ps & nt(tst) >= xo;
    N = struct('norm', nrm, 'baseline', b, 'threshold', thr, 'crossOct', xo, 'crossQ', opts.CrossQ, ...
        'falseTestS', sum(ps)/(numel(tst)*W), 'falseTestW', mean(ps), ...
        'crossTestS', sum(px)/(numel(tst)*W), 'crossTestW', mean(px), 'bandOcc', occ/max(sum(occ), 1), ...
        'alpha', opts.Alpha, 'windowS', W, 'calScores', sc(cal), 'calNet', cn, ...
        'testScores', sc(tst), 'testSpan', sp(tst), 'testSweep', sw(tst), 'testNet', nt(tst));
end

% Author: Diellor Basha, 2026
