function N = jointpathnull(Ysur, G, nF, rate, opts)
% DETECT.JOINTPATHNULL  Normalisation, baseline, threshold and held-out false rates for rheome.detect.jointpath.
%
%   N = rheome.detect.jointpathnull(Ysur, G, nF, rate, LogF=log2fc)              % Alpha 0.1/s, BaseQ 0.9
%   N = rheome.detect.jointpathnull(Ysur, G, nF, rate, LogF=log2fc, Alpha=0.05, CrossOct=1)
%
% Ysur is a NULL band-envelope field -- a surrogate of the data that keeps each (band, place) series'
% spectrum and destroys its alignment with the others -- cut into windows of nF frames and concatenated:
% [nTile x nBand x nF*nWin]. As rheome.detect.tilepathnull: odd windows CALIBRATE, even windows TEST.
%   norm       [nTile x nBand] median of each series over the calibration windows (jointpath's Norm)
%   baseline   the BaseQ quantile of the normalised calibration values
%   threshold  the best-path score at most Alpha*WindowS calibration windows exceed (Alpha false paths/s)
%   falseTestS passing best paths per second on the TEST windows (any path)
%   crossTestS the same, counting only paths whose frequency centroid spans >= CrossOct octaves: the false
%              rate of the claim this method exists to make, "a pattern moved across frequency"
%
% OUTPUT (struct N): .norm .baseline .threshold .falseTestS .crossTestS .alpha .windowS .crossOct
%   .calScores .testScores .testSpan .testSweep [nWin/2 x 1]
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
        opts.CrossOct (1,1) double {mustBePositive} = 1
    end
    [nTile, nB, nT] = size(Ysur);  nWin = floor(nT/nF);
    if nWin < 4, error('detect:jointpathnull:short', 'Need >= 4 null windows, got %d.', nWin); end
    W = opts.WindowS;  if isempty(W), W = nF/rate; end
    cal = 1:2:nWin;  tst = 2:2:nWin;
    ci = reshape((cal(:) - 1)*nF + (1:nF), 1, []);
    nrm = median(Ysur(:, :, ci), 3);  nrm(~(nrm > 0)) = realmin;      % a dead series must not divide by 0
    Zc = Ysur(:, :, ci) ./ nrm;  b = quantile(Zc(:), opts.BaseQ);
    sc = zeros(nWin, 1);  sp = nan(nWin, 1);  sw = sp;
    for w = 1:nWin
        p = rheome.detect.jointpath(Ysur(:, :, (w-1)*nF + (1:nF)), G, Norm=nrm, Baseline=b, LogF=opts.LogF, SampleRate=rate);
        if p.nTracks, sc(w) = p.tracks(1).score;  sp(w) = p.tracks(1).octSpan;  sw(w) = p.tracks(1).sweepOctS; end
    end
    s = sort(sc(cal), 'descend');  k = floor(opts.Alpha*W*numel(s));
    if k >= numel(s), thr = 0; else, thr = s(k + 1); end
    ps = sc(tst) > thr;
    N = struct('norm', nrm, 'baseline', b, 'threshold', thr, 'falseTestS', sum(ps)/(numel(tst)*W), ...
        'crossTestS', sum(ps & sp(tst) >= opts.CrossOct)/(numel(tst)*W), 'alpha', opts.Alpha, 'windowS', W, ...
        'crossOct', opts.CrossOct, 'calScores', sc(cal), 'testScores', sc(tst), 'testSpan', sp(tst), 'testSweep', sw(tst));
end

% Author: Diellor Basha, 2026
