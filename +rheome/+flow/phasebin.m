function out = phasebin(A, phase, nBins)
% FLOW.PHASEBIN  Collapse a coefficient series onto a common PHASE axis.
%
%   out = rheome.flow.phasebin(A, phase, nBins)
%
% ⭐ THIS IS WHAT MAKES MULTIRATE RECTANGULAR AGAIN. Under a multirate scheme every filter
% runs at its own sample rate (r_m = samplesPerCycle * f_m), so the TIME axes disagree and
% the coefficients cannot be stacked into one array. Binning by the reference rhythm's
% PHASE gives every filter the same nBins per cycle -- and the results stack:
%
%     A_m [Ks x nT_m]   (nT_m differs per filter)
%          -> phasebin ->   [Ks x nBins]   (identical for every filter)
%
% At 30 bins, 800 modes and 7 filters that is 1.3 MB against a 1.16 GB full-rate tensor.
%
% ⚠ THE TWO ACCUMULATORS ANSWER DIFFERENT QUESTIONS. Do not substitute one for the other:
%
%   .mean    COMPLEX mean per bin. Keeps only what is PHASE-LOCKED to the reference --
%            the cycle-averaged flow map, the thing you animate. Activity not locked to
%            the reference CANCELS, by design.
%   .energy  mean |A|^2 per bin. Keeps ALL power at that phase, locked or not; nothing
%            can cancel. Use it to ask WHEN in the cycle the flow is strongest.
%
% Their ratio, |mean|^2 / energy, is a phase-locking measure in [0,1].
%
% ⚠ EMPTY BINS ARE NaN, NOT ZERO. A bin no sample reached is UNKNOWN. Returning zero would
% let "no measurement" average into downstream statistics as "no flow" -- silently, and in
% the direction of weakening any real effect.
%
% ⚠ PHASE IS TAKEN ON (-pi, pi] AND WRAPS. -pi and +pi are the SAME phase and land in the
% same bin; treating them as distinct creates a phantom bin at the wrap point.
%
% INPUTS:
%   A       [nMode x nT] coefficients (complex; real is accepted)
%   phase   [1 x nT] reference phase per sample, radians
%   nBins   number of phase bins (the multirate samplesPerCycle is the natural choice)
%
% OUTPUT (struct out):
%   .mean    [nMode x nBins] complex   .energy [nMode x nBins] real
%   .count   [1 x nBins]               .phase  [1 x nBins] bin-centre phase (rad)
%   .locking [nMode x nBins] |mean|^2 ./ energy, in [0,1]
%
% See also: rheome.flow.rateplan, rheome.flow.demod, rheome.flow.observables
%
% Author: Diellor Basha, 2026

    nT = size(A, 2);
    if numel(phase) ~= nT
        error('flow:phasebin:size', ...
            'phase has %d samples but A has %d columns.', numel(phase), nT);
    end
    if ~isscalar(nBins) || nBins < 1 || mod(nBins,1) ~= 0
        error('flow:phasebin:bins', 'nBins must be a positive integer.');
    end

    nMode = size(A, 1);
    phase = double(phase(:)).';

    % Wrap to (-pi, pi], then map to 1..nBins. mod() puts +pi and -pi on the same point,
    % which is what makes the wrap seam a single bin rather than two half-bins.
    w   = mod(phase + pi, 2*pi);                     % [0, 2*pi)
    idx = floor(w / (2*pi) * nBins) + 1;
    idx = min(max(idx, 1), nBins);                   % guard the closed upper edge

    cnt = accumarray(idx(:), 1, [nBins 1]).';

    sumC = complex(zeros(nMode, nBins));
    sumE = zeros(nMode, nBins);
    for b = 1:nBins
        m = (idx == b);
        if ~any(m), continue; end
        Ab = A(:, m);
        sumC(:, b) = sum(Ab, 2);
        sumE(:, b) = sum(abs(Ab).^2, 2);
    end

    denom = cnt;  denom(denom == 0) = NaN;           % empty bin -> NaN, never 0
    out.mean   = sumC ./ denom;
    out.energy = sumE ./ denom;
    out.count  = cnt;
    out.phase  = -pi + ((1:nBins) - 0.5) * (2*pi/nBins);
    out.locking = abs(out.mean).^2 ./ out.energy;
end

% Author: Diellor Basha, 2026
