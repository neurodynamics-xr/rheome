function out = vortexspectrum(name, varargin)
% FLOW.VORTEXSPECTRUM  Vortices in several bands on a 1/f background, matched to a real sensor PSD.
%
%   out = rheome.flow.vortexspectrum('sub01')
%   out = rheome.flow.vortexspectrum(name, Bands=[2 4; 8 16; 16 32], Duration=20)
%
% ⭐⭐ THE AMPLITUDE IS SOLVED, NOT GUESSED. The forward is linear, so a band's sensor power scales
% as the square of the source moment. Measuring the band's power at unit moment and dividing once
% gives the moment that reproduces the measured PERIODIC power exactly -- no sweep, no fitting. That
% is the whole answer to "how do we get a spectral signature matching the data": do not set source
% amplitudes and hope, invert the forward on the band powers.
%
% THE CONSTRUCTION, in two pieces that add at the sensors:
%   background  synthesised IN the Dirac mode basis at the fitted exponent (rheome.spectral.synthesise), so
%               it is in the span by construction, and forwarded as Gmode*c_ap
%   rhythms     one TRAIN of rheome.flow.vortexatom per band, at the tile stride, each with a random phase
% ⭐ A single atom cannot carry a band's power: it is time-localised, and the real rhythm runs the
% whole record. A train at strideSelect = supportSec is the tiling that fills the band -- the same
% lattice rheome.selection.wavelettile navigates -- so the band is continuously occupied by rotating
% sources rather than by one burst.
%
% ⭐⭐ THE BACKGROUND MATCHES; THE RHYTHMS ARE WHAT COST ACCURACY. Split by region, the
% channel-averaged log-PSD correlation is +0.96 OUTSIDE the banded range, where only the background
% acts, and +0.63 INSIDE 2-32 Hz where the atoms live (0.110 against 0.241 decades). A constant-Q
% atom delivers its band's power as one octave-wide bump at fc, which is not the data's within-band
% shape, so adding rhythms to a good background makes the curve worse while making the BAND POWERS
% right: 1.021 / 1.020 / 1.020 / 0.994 over delta / theta / alpha / beta. ⚠ Decide which you need.
% A background alone reproduces the spectrum better and attributes nothing; this attributes each
% band's power to a rotating source and pays for it in spectral detail.
%
% ⚠ WHAT IS MATCHED AND WHAT IS NOT. The channel-averaged PSD shape and the per-band power are
% matched. The SPATIAL distribution is not: the background is synthesised with equal power per mode,
% which is white across the Dirac spectrum and is not the brain's spatial covariance, and each band's
% rhythm sits at ONE location. So this reproduces a spectral signature, not a full recording. ⭐ It
% is the right object for asking whether a band's measured power is consistent with a rotating
% source of a given scale; it is the wrong object for spatial-covariance work.
%
% ⚠ LOW BANDS NEED LONG RECORDS. An atom's support is 15.05/fc, so delta at 2.83 Hz supports 5.3 s
% and needs several times that. Bands whose support exceeds a third of Duration are skipped with a
% warning rather than silently truncated.
%
% INPUTS
%   Bands      [nB x 2] Hz, default [2 4; 4 8; 8 16; 16 32]
%   Duration   seconds (20)             FitRange  [1 45] for the aperiodic fit
%   WavelengthMM (140)                  Vertices  one seed per band, or [] to spread them
%   Hemi ("L")  Bases  Gauge  Study  Dirac  Verbose
%
% OUTPUT (struct out)
%   .B        [nCh x nT] the synthesised sensor record       .Breal the real one, same window
%   .f .Psim .Preal   channel-averaged Welch spectra
%   .bands    table: band, fLo, fHi, fc, periodicFrac, momentNAm, nAtoms, supportSec, used
%   .chi .knee  the fitted aperiodic parameters of the REAL data
%   .check    logCorr, medianAbsLogRatio, chiReal, chiSim, bandPowerRatio
%
% See also: rheome.spectral.aperiodic, rheome.spectral.synthesise, rheome.flow.vortexatom, rheome.forward.diracgain,
%           rheome.selection.wavelettile
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Bands', [2 4; 4 8; 8 16; 16 32], @(x) size(x,2)==2);
    p.addParameter('Duration', 20, @isscalar);
    p.addParameter('FitRange', [1 45], @(x) numel(x)==2);
    p.addParameter('WavelengthMM', 140, @isscalar);
    p.addParameter('Vertices', [], @isnumeric);
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Gauge', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Study', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Dirac', [], @(x) isempty(x) || isstruct(x));
    % ⭐ "asfitted" is right ONCE FitIter/FitThresh make the fit a floor. "floor" rescales the
    %   background to stay under the observed PSD and was needed only with rheome.spectral.aperiodic's
    %   defaults; it is driven by a single worst bin, weakens the background everywhere, and
    %   measured it makes the match worse (r +0.79 against +0.83, 0.362 decades against 0.196).
    p.addParameter('BackgroundFit', "asfitted", @(x) any(strcmpi(x, ["floor","asfitted"])));
    % ⭐ NOT rheome.spectral.aperiodic's defaults. At nIter 3 / thresh 2 the fit sits at the MIDDLE of the
    %   data (49.8% of bins above it) and cannot serve as a background. At 6 / 0.5 it is over at
    %   only 5.3% of bins with R2 0.992 and chi 2.40. ⚠ It is not monotone in either knob: 12 / 0.5
    %   collapses to chi 0.80 with 50.8% over, so these are a measured sweet spot, not "more is
    %   better", and changing them means re-measuring the overshoot.
    p.addParameter('FitIter', 6, @isscalar);
    p.addParameter('FitThresh', 0.5, @isscalar);
    p.addParameter('Verbose', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    B = o.Bases;  if isempty(B), B = rheome.load.bases(name); end
    Hm = B.(char(o.Hemi));  S = Hm.S;  nVh = size(S.Vertices,1);  gv = double(Hm.gv(:))';
    g = o.Gauge;  if isempty(g), g = rheome.operators.gauge(S.Vertices, double(S.Faces), Method="diffusion"); end
    st = o.Study; if isempty(st), st = rheome.load.study(name); end
    dd = o.Dirac; if isempty(dd), dd = rheome.load.dirac(name); end
    [G, idx] = rheome.forward.leadfield(st, GlobalVertices=gv);

    %% 0. the real record, decimated to the working rate
    fsr = st.rec.sfreq;  dec = max(1, round(fsr/300));  fs = fsr/dec;
    nT  = round(o.Duration*fs);
    Fr  = double(st.rec.F(idx, 1:dec:end));
    assert(size(Fr,2) >= nT, 'rheome.flow.vortexspectrum: the record is shorter than Duration.');
    Breal = Fr(:, 1:nT);

    %% 1. fit the REAL sensor background
    win = hann(round(4*fs));
    [Pr, f] = pwelch(Breal.', win, [], [], fs);
    Pr = mean(Pr, 2);                                    % channel-averaged
    kf = f > o.FitRange(1) & f < o.FitRange(2);
    ap = rheome.spectral.aperiodic(Pr(kf), f(kf), ...
             struct('knee', true, 'nIter', o.FitIter, 'thresh', o.FitThresh));
    df = mean(diff(f));

    %% 2. per-band periodic power, from the fit
    nB = size(o.Bands,1);
    Pap = nan(size(f));  Pap(kf) = ap.Pap;
    per = zeros(nB,1);  tgt = zeros(nB,1);  fcb = zeros(nB,1);
    for b = 1:nB
        kb = f >= o.Bands(b,1) & f <= o.Bands(b,2) & kf;
        if ~any(kb), continue; end
        ex = max(Pr(kb) - Pap(kb), 0);
        per(b) = sum(ex) / max(sum(Pr(kb)), realmin);
        tgt(b) = sum(ex) * df;                            % provisional; corrected in step 4
        fcb(b) = sqrt(o.Bands(b,1)*o.Bands(b,2));
    end

    %% 3. the 1/f background, synthesised in the Dirac mode basis
    cols = find(dd.Hemisphere == 1);
    Phi  = dd.Phi(reshape((gv'-1)*4 + (1:4), [], 1), cols);
    Gm   = rheome.forward.diracgain(G, struct('Phi', Phi, 'nVert', nVh, 'nModes', numel(cols)));
    % ⚠⚠ USE THE FITTED CURVE, NOT ITS PARAMETERS. rheome.spectral.aperiodic searches the knee on a grid
    %   topping out at 1000, and this data pins it there, so a background built from (chi, knee) is
    %   not free to match the low end -- measured, it came out 20% hot in delta and left no room for
    %   a delta rhythm at all. Passing ap.Pap itself reproduces the background that was fitted.
    Cap  = rheome.spectral.synthesise(nT, fs, NaN, ...
              struct('nS', numel(cols), 'shape', ap.Pap, 'shapeF', f(kf)));
    Bbg  = Gm*Cap;
    % one global scale so the background carries the measured APERIODIC power
    Pbg  = mean(pwelch(Bbg.', win, [], [], fs), 2);
    scB  = sqrt(sum(Pap(kf)) / max(sum(Pbg(kf)), realmin));
    % ⚠⚠ THE FITTED BACKGROUND IS NOT A FLOOR, and a synthesis needs one. Measured on this
    %   recording, ap.Pap EXCEEDS the observed PSD at 49.8% of in-fit bins and by band gives
    %   Pap/Preal of 1.101 / 1.034 / 0.345 / 0.745 / 1.020 over 2-4 / 4-8 / 8-16 / 16-32 / 32-45 Hz.
    %   specparam's peak rejection down-weights points ABOVE the fit but still lands near the
    %   middle of the data, so "observed minus background" is NEGATIVE in delta, theta and gamma
    %   and those bands get no rhythm at all. "floor" rescales the background to the largest
    %   multiple of its own shape that stays under the observed PSD, which is the definition a
    %   synthesis actually needs. ⚠ It makes the background weaker than the fit, so the periodic
    %   fractions reported alongside it are correspondingly larger -- a modelling choice, not a
    %   measurement of how much of the brain's power is rhythmic.
    if strcmpi(o.BackgroundFit, "floor")
        Pbg1 = mean(pwelch((Bbg*scB).', win, [], [], fs), 2);
        scB  = scB * sqrt(min(Pr(kf) ./ max(Pbg1(kf), realmin)));
    end
    Bbg  = Bbg*scB;
    % ⭐ TOP UP AGAINST THE BACKGROUND THAT WAS ACTUALLY BUILT, not the fitted curve. One global
    %   scale makes the background match Pap in TOTAL over the fit range but not band by band, and
    %   sizing each rhythm against the fit instead left every band 2-21% hot. Measuring what the
    %   synthesised background delivers in each band and asking the rhythm for the remainder makes
    %   the band powers exact by construction.
    Pbg2 = mean(pwelch(Bbg.', win, [], [], fs), 2);
    for b = 1:nB
        kb = f >= o.Bands(b,1) & f <= o.Bands(b,2) & kf;
        if ~any(kb), continue; end
        tgt(b) = max(sum(Pr(kb)) - sum(Pbg2(kb)), 0) * df;
    end

    %% 4. one train of vortex atoms per band, moment SOLVED from the target power
    if isempty(o.Vertices)
        vseed = round(linspace(0.15, 0.85, nB)*nVh);
    else
        vseed = o.Vertices(:)';
    end
    Bsig = zeros(size(G,1), nT);
    mom  = nan(nB,1);  nAt = zeros(nB,1);  sup = nan(nB,1);  used = false(nB,1);
    rng(19);
    for b = 1:nB
        at = rheome.flow.vortexatom(name, Vertex=vseed(b), WavelengthMM=o.WavelengthMM, ...
                 Band=o.Bands(b,:), Duration=o.Duration, SampleRate=fs, ...
                 Hemi=o.Hemi, Bases=B, Gauge=g, Check=false);
        sup(b) = at.check.supportSec;
        if sup(b) > o.Duration/3
            warning('flow:vortexspectrum:shortRecord', ...
                ['band %g-%g Hz supports %.2f s in a %.1f s record; skipping. Lengthen ' ...
                 'Duration to at least %.0f s for this band.'], o.Bands(b,1), o.Bands(b,2), ...
                sup(b), o.Duration, 3*sup(b));
            continue;
        end
        used(b) = true;
        pA = G*at.JA;  pB = G*at.JB;
        % ⭐ the train: atoms every supportSec, each with its own phase, filling the band
        stride = round(sup(b)*fs);
        starts = 1:stride:nT;
        nAt(b) = numel(starts);
        aT = zeros(1,nT);  bT = zeros(1,nT);
        for j = 1:numel(starts)
            sh  = starts(j) - 1;
            psh = circshift(at.psi, sh) * exp(1i*2*pi*rand);
            aT  = aT + real(psh);  bT = bT + imag(psh);
        end
        U  = pA*aT + pB*bT;                               % the band at UNIT moment
        Pu = mean(pwelch(U.', win, [], [], fs), 2);
        kb = f >= o.Bands(b,1) & f <= o.Bands(b,2);
        sc = sqrt(tgt(b) / max(sum(Pu(kb))*df, realmin));  % ⭐ one division, exact
        Bsig = Bsig + U*sc;
        mom(b) = sc * sum(vecnorm(reshape(at.JA,3,[])',2,2)) * 1e9;
        if o.Verbose
            fprintf('  %5.1f-%-5.1f Hz  periodic %.3f  %2d atoms  moment %.2f nAm\n', ...
                o.Bands(b,1), o.Bands(b,2), per(b), nAt(b), mom(b));
        end
    end

    %% 5. the synthesised record, and how well it matches
    Bs = Bbg + Bsig;
    Ps = mean(pwelch(Bs.', win, [], [], fs), 2);
    kc = f >= o.FitRange(1) & f <= o.FitRange(2);
    apS = rheome.spectral.aperiodic(Ps(kc), f(kc), ...
              struct('knee', true, 'nIter', o.FitIter, 'thresh', o.FitThresh));

    out.B = Bs;  out.Breal = Breal;  out.f = f;  out.Psim = Ps;  out.Preal = Pr;
    out.chi = ap.exponent;  out.knee = ap.knee;  out.fs = fs;  out.nT = nT;
    out.bands = table((1:nB)', o.Bands(:,1), o.Bands(:,2), fcb, per, mom, nAt, sup, used, ...
        'VariableNames', {'band','fLo','fHi','fc','periodicFrac','momentNAm','nAtoms','supportSec','used'});
    out.check.logCorr           = corr(log(Ps(kc)), log(Pr(kc)));
    out.check.medianAbsLogRatio = median(abs(log10(Ps(kc)./Pr(kc))));
    out.check.chiReal           = ap.exponent;
    out.check.chiSim            = apS.exponent;
    bp = nan(nB,1);
    for b = 1:nB
        if ~used(b), continue; end
        kb = f >= o.Bands(b,1) & f <= o.Bands(b,2);
        bp(b) = sum(Ps(kb)) / max(sum(Pr(kb)), realmin);
    end
    out.check.bandPowerRatio = bp;
    out.check.fitOverFrac    = mean(ap.Pap > Pr(kf));   % ⚠ a background should be near zero here
    out.check.fitR2          = ap.r2;
end

% Author: Diellor Basha, 2026
