function N = tilepathnull(Ysur, G, nF, rate, opts)
% DETECT.TILEPATHNULL  Baseline and score threshold for rheome.detect.tilepath from null windows.
%
%   N = rheome.detect.tilepathnull(Ysur, G, nF, rate)                       % Alpha 0.1/s, BaseQ 0.9
%   N = rheome.detect.tilepathnull(Ysur, G, nF, rate, Alpha=0.05, WindowS=2)
%
% Ysur is a NULL tile field -- noise alone, or a surrogate of the data -- cut into observation
% windows of nF frames and concatenated: [nTile x nF*nWin]. Odd windows CALIBRATE and even windows
% TEST, so the false rate reported is measured on windows the thresholds never saw:
%   baseline   the BaseQ quantile of the calibration windows' tile values: what a walk must beat to
%              gain score (rheome.detect.tilepath's Baseline)
%   threshold  the best-path score that at most Alpha*WindowS calibration windows exceed -- a false
%              rate of Alpha paths per second
%
% ⭐ A RATE PER SECOND, NOT A PER-WINDOW p-VALUE: a reader of a trajectory map needs to know how many
% of its paths per second are noise. One best path per window is scored.
%
% OUTPUT (struct N)
%   .baseline .threshold .falseTestS (measured on the test windows) .alpha .windowS
%   .calScores .testScores [nWin/2 x 1]
%
% See also: rheome.detect.tilepath, rheome.detect.tracknull
%
% Author: Diellor Basha, 2026

    arguments
        Ysur double
        G struct
        nF (1,1) double {mustBeInteger, mustBePositive}
        rate (1,1) double {mustBePositive}
        opts.Alpha (1,1) double {mustBePositive} = 0.1
        opts.BaseQ (1,1) double {mustBeInRange(opts.BaseQ, 0, 1)} = 0.9
        opts.WindowS = []
    end
    nWin = floor(size(Ysur, 2) / nF);
    if nWin < 4, error('detect:tilepathnull:short', 'Need >= 4 null windows, got %d.', nWin); end
    W = opts.WindowS;  if isempty(W), W = nF / rate; end
    cal = 1:2:nWin;  tst = 2:2:nWin;
    calv = cell2mat(arrayfun(@(w) reshape(Ysur(:, (w-1)*nF + (1:nF)), [], 1), cal(:), 'uni', 0));
    b = quantile(calv, opts.BaseQ);
    sc = zeros(nWin, 1);
    for w = 1:nWin
        p = rheome.detect.tilepath(Ysur(:, (w-1)*nF + (1:nF)), G, Baseline=b, SampleRate=rate);
        if p.nTracks, sc(w) = p.tracks(1).score; end
    end
    s = sort(sc(cal), 'descend');  k = floor(opts.Alpha * W * numel(s));
    if k >= numel(s), thr = 0; else, thr = s(k + 1); end
    N = struct('baseline', b, 'threshold', thr, 'falseTestS', sum(sc(tst) > thr) / (numel(tst) * W), ...
        'alpha', opts.Alpha, 'windowS', W, 'calScores', sc(cal), 'testScores', sc(tst));
end

% Author: Diellor Basha, 2026
