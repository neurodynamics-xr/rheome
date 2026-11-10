function S = burststages(z, good, fps, per, opts)
% CORES.BURSTSTAGES  The frames of each burst stage: onset (build-up), within (plateau), offset (decay), and the fits.
%
%   S = rheome.cores.burststages(z, good, fps, per)
%
% z [C x m] complex band frames (rheome.cores.bandframes), good [1 x m] logical, fps frames a second, per frames per cycle.
% Pipeline (nxr's constant-Q form of the lab's standard): amplitude a = sqrt(mean over channels |z|^2), smoothed 1 cycle
% (rheome.cores.boxsmooth, zero-padded) -> candidate bursts (rheome.cores.burstenvelope) -> each located by rheome.cores.burstfit
% (t0 = start/fps, t1 = stop/fps; margin 10 cycles, smoothing 1 cycle, tau <= 10 cycles; a cycle = per/fps s). A fit is kept when ok,
% not at its bound and r2 >= MinR2 (0.7). The onset mask is [tOn +- 2 tauR], within [tOn + 2 tauR, tOff - 2 tauF], offset
% [tOff +- 2 tauF]: frames lo = max(0, ceil(iv0*fps)) ... min(m, floor(iv1*fps) + 1) (0-based, END INCLUSIVE of the last frame <= iv1),
% ANDed with good last. A short transient may pass the 3-cycle gate (the smoothings widen it) and is fitted: that is nxr's behaviour.
%
% OUTPUT S (struct): .onset .within .offset [1 x m] logical; .fits [nkept x 5] = tOn tOff tauR tauF r2 of the kept fits;
% .candidates [n x 2] 1-based inclusive frames (burstenvelope); .all [n x 1] struct array of every candidate's burstfit result
% (kept or not); .amp1 the 1-cycle smoothed amplitude; .thresholds [median, robust SD, hi, lo].
%
% Ports nxr-analysis src/nxr/vortex.py burst_stages (line 845, commit 30b9622).
%
% See also: rheome.cores.bandframes, rheome.cores.burstenvelope, rheome.cores.burstfit, rheome.cores.burstmask
%
% Author: Diellor Basha, 2026
    arguments
        z (:,:) {mustBeNumeric}
        good (1,:) logical
        fps (1,1) double {mustBePositive}
        per (1,1) double {mustBeInteger, mustBePositive}
        opts.MinR2 (1,1) double = 0.7
    end
    m = size(z, 2);  cyc = per / fps;
    a = sqrt(mean(abs(double(z)).^2, 1));
    a = rheome.cores.boxsmooth(a, per);
    [cand, th] = rheome.cores.burstenvelope(a, good, per);
    names = {'onset', 'within', 'offset'};
    for q = 1:3, S.(names{q}) = false(1, m); end
    S.fits = zeros(0, 5);  allfits = cell(0, 1);
    for c = 1:size(cand, 1)
        i = cand(c, 1) - 1;  j = cand(c, 2);                       % 0-based [start, stop)
        R = rheome.cores.burstfit(a, fps, i / fps, j / fps, Margin=10*cyc, Smooth=cyc, TauMax=10*cyc);
        allfits{c, 1} = R;
        if ~R.ok || R.atBound || R.r2 < opts.MinR2, continue; end
        S.fits(end+1, :) = [R.tOn R.tOff R.tauR R.tauF R.r2]; %#ok<AGROW>
        ivs = {R.buildup, R.plateau, R.decay};
        for q = 1:3
            iv = ivs{q};
            if isempty(iv), continue; end
            lo = max(0, ceil(iv(1) * fps));  hi = min(m, floor(iv(2) * fps) + 1);
            if hi > lo, S.(names{q})(lo+1:hi) = true; end
        end
    end
    for q = 1:3, S.(names{q}) = S.(names{q}) & good; end
    S.candidates = cand;  S.amp1 = a;  S.thresholds = th;
    if isempty(allfits), S.all = struct('tOn', {}); else, S.all = vertcat(allfits{:}); end
end
