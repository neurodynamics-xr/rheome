function out = phaselock(z, refPhase, nBins, varargin)
% FLOW.PHASELOCK  Is the cortical field reproducibly organised by a reference rhythm's phase?
%
%   out = rheome.flow.phaselock(z, refPhase, nBins)
%
% Pools a complex cortical field over many cycles of an external reference -- the sensor-space
% alpha phase -- and asks whether the cycle-averaged spatial map is REAL, in the only sense
% that can be checked from data: whether it reproduces on cycles it was not built from.
%
% ⭐ THE OPERATOR IS COMMON-MODE HERE, WHICH IS WHY THIS WORKS WHERE A CORE CENSUS DOES NOT.
% The inverse smooths every frame identically, so it cannot manufacture a DIFFERENCE between
% one reference phase and another. Anything that varies with reference phase survived the
% operator rather than being made by it. Absolute spatial quantities -- wavelength, core
% count, propagation speed -- do NOT have this property: measured on real alpha they matched
% a phase-randomised surrogate to within a few percent, because the operator sets them.
%
% ⚠ PLV IS BIASED UPWARD BY FINITE SAMPLE SIZE, so it must not be the headline. With C cycles
% and no locking whatsoever, E|PLV| ~ 1/sqrt(C): at 20 cycles that is 0.22 of apparent
% locking from nothing at all. .splitHalf is built from disjoint halves of the cycles and has
% no such bias -- it is ~0 under the null however few cycles there are. Report it, and use
% .plv only for the spatial pattern once .splitHalf says there is one.
%
% ⚠ STANDING vs TRAVELLING IS NOT A RANK QUESTION, and reading it off the binned map's
% singular values is simply wrong. At one frequency z(v,t) = w(v)*exp(i*omega*t), so the
% binned map is ONE complex spatial pattern times ONE phase function -- rank 1 BY
% CONSTRUCTION, however the pattern travels. Measured: a deliberately planted travelling
% rotor drove the binned map's "travel index" DOWN, 0.114 -> 0.023, as its amplitude rose.
% The travelling structure lives in the spatially-varying PHASE of w, so .standing is the
% eigenvalue split of the [nV x 2] real matrix [Re w, Im w]: all vertices on one line
% through the origin means one common phase and a standing pattern.
%
% ⚠ CYCLES, NOT SAMPLES, ARE THE UNIT OF EVIDENCE. Neighbouring samples within one cycle are
% not independent, so splitting by sample would put near-duplicates on both sides and inflate
% the reproducibility toward 1. Splitting alternates whole CYCLES.
%
% INPUTS:
%   z         [nV x nT] complex cortical field (the analytic curl)
%   refPhase  [1 x nT]  reference phase per sample, radians, on (-pi, pi]
%   nBins     phase bins per cycle (30 is the project's standard)
%   'MinCycle'  cycles shorter than this fraction of the median are dropped as phase-slips.
%               Default 0.5.
%
% OUTPUT (struct out):
%   .splitHalf  [nV x 1]  complex coherence between odd- and even-cycle binned maps, per
%                         vertex, real part. ~0 under the null at ANY number of cycles.
%   .global     scalar    the same coherence taken over all vertices at once
%   .plv        [nV x 1]  |mean|/sqrt(energy) pooled over all cycles (BIASED -- see above)
%   .pref       [nV x 1]  preferred reference phase, radians
%   .binned     [nV x nBins] complex cycle-averaged map
%   .locked     [nV x 1]  the phase-locked COMPLEX pattern (first Fourier mode over bins)
%   .standing   scalar in [0.5, 1]: 1 = every vertex peaks at the same reference phase
%                         (a standing pattern), 0.5 = phases spread over the full plane
%                         (a travelling / rotating one)
%   .nCycle .cycleStart
%
% See also: rheome.flow.phasebin, rheome.flow.phasegradient, rheome.detect.chargedensity
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('MinCycle', 0.5);
    p.parse(varargin{:});
    o = p.Results;

    refPhase = refPhase(:).';
    if size(z,2) ~= numel(refPhase)
        error('flow:phaselock:size', ...
            'z has %d samples but refPhase has %d.', size(z,2), numel(refPhase));
    end
    if isreal(z)
        error('flow:phaselock:real', ...
            'z is real -- this needs the ANALYTIC field. Apply hilbert() in time first.');
    end

    % Cycle boundaries: the reference phase wrapping from + back to -.
    d  = diff(refPhase);
    st = [1, find(d < -pi) + 1, numel(refPhase)+1];
    len = diff(st);
    keep = len >= o.MinCycle * median(len);           % drop phase-slip fragments
    st  = st(1:end-1);
    st  = st(keep);  len = len(keep);
    nC  = numel(st);
    if nC < 4
        error('flow:phaselock:cycles', ...
            'Only %d usable cycles; the split-half statistic needs at least 4.', nC);
    end

    half = false(1, numel(refPhase));                  % alternate whole CYCLES, not samples
    for k = 1:2:nC
        half(st(k) : st(k)+len(k)-1) = true;
    end

    A = rheome.flow.phasebin(z,          refPhase,          nBins);
    B = rheome.flow.phasebin(z(:,half),  refPhase(half),    nBins);
    C = rheome.flow.phasebin(z(:,~half), refPhase(~half),   nBins);

    out.binned = A.mean;
    out.plv    = sqrt(sum(abs(A.mean).^2, 2, 'omitnan') ./ ...
                      max(sum(A.energy, 2, 'omitnan'), realmin));
    out.pref   = angle(sum(A.mean, 2, 'omitnan'));

    % Coherence between the two halves' binned maps: <B, C> / sqrt(<B,B><C,C>). Real part,
    % so a map that reproduces with a PHASE SHIFT does not count as reproducing.
    ok = isfinite(B.mean) & isfinite(C.mean);
    Bm = B.mean;  Cm = C.mean;  Bm(~ok) = 0;  Cm(~ok) = 0;
    num = sum(Bm .* conj(Cm), 2);
    den = sqrt(sum(abs(Bm).^2, 2) .* sum(abs(Cm).^2, 2));
    out.splitHalf = real(num ./ max(den, realmin));

    gn = sum(Bm(:) .* conj(Cm(:)));
    gd = sqrt(sum(abs(Bm(:)).^2) * sum(abs(Cm(:)).^2));
    out.global = real(gn / max(gd, realmin));

    % The phase-locked complex pattern: the first Fourier mode across phase bins.
    th = A.phase(:).';
    Am = A.mean;  Am(~isfinite(Am)) = 0;
    w  = Am * exp(-1i*th).' / nBins;
    out.locked = w;
    W  = [real(w), imag(w)];
    sv = svd(W, 'econ');
    out.standing = sv(1)^2 / max(sv(1)^2 + sv(2)^2, realmin);

    out.nCycle     = nC;
    out.cycleStart = st;
end

% Author: Diellor Basha, 2026
