function d = dispersion(obj, C)
% DISPERSION  Fit the speed of a travelling structure from the joint spectrum.
%
%   d = dispersion(jfb, C)
%
% ⭐ WHY THIS IS THE USEFUL READOUT. A structure travelling at speed c satisfies
% omega = c*sqrt(lambda), so a constant-speed feature appears as a DIAGONAL RIDGE in the
% (sqrt(lambda), omega) plane. Reading its slope MEASURES the speed, which turns a
% speed-selective filter from an assumed velocity into one fitted to the data.
%
% ⭐ THE MEMBERS ARE NOT INVOLVED. The ridge is taken from the FULL (lambda, omega)
% energy: with a handful of overlapping members the wavenumber axis is far too coarse to
% fit a slope, and doing so returns nonsense.
%
% ⭐ THE FIT IS ENERGY-WEIGHTED. Unweighted, every bin counts equally, so bins carrying
% no signal -- outside the diagonal's span, where kbar is merely the centroid of the
% lambda axis -- dominate by sheer number and bias the slope badly.
%
% ⚠ THE BAND MUST SPAN THE DIAGONAL. The slope is recoverable only if the retained
% frequency axis covers c*sqrt(lambda)/2pi across the wavenumbers present. A narrow band
% cannot measure speed this way: an 8-13 Hz retention sees a 5 Hz slice of a diagonal
% spanning tens of Hz. CHECK .fSupport against the diagonal you expect before believing
% .slope.
%
% OUTPUT (struct d):
%   .slope    fitted speed (operator length units per second)
%   .slopeR2  weighted coefficient of determination
%   .kbar     [1 x nOk] energy-weighted mean wavenumber per frequency -- the ridge itself
%   .fRidge   [1 x nOk] the frequencies it was evaluated at
%   .wRidge   [1 x nOk] the normalised weights used
%   .fSupport [f5 f95] the frequency interval carrying 90% of the energy
%
% See also: scalogram, rheome.jointfilterbank.speedkernel
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda;  f = obj.Frequencies;
    K = numel(lam);  nO = numel(f);
    if ~isequal(size(C), [K nO])
        error('jointfilterbank:size', ...
            'C is %s but the bank grid is %s.', mat2str(size(C)), mat2str([K nO]));
    end

    P  = abs(C).^2;
    kk = sqrt(lam);
    w  = sum(P, 1);
    ok = w > 0;

    d = struct('slope', NaN, 'slopeR2', NaN, 'kbar', [], 'fRidge', [], ...
               'wRidge', [], 'fSupport', [NaN NaN]);
    if nnz(ok) <= 2, return; end

    kbar = (kk.' * P(:,ok)) ./ w(ok);          % [1 x nOk] mean wavenumber per frequency
    om   = 2*pi*f(ok);

    % ⚠ NO SPREAD IN kbar MEANS NO DIAGONAL TO FIT. This happens when the ridge lies
    % almost entirely off the retained axis, so every frequency reports the same mean
    % wavenumber. Return NaN rather than let the least-squares go rank-deficient and
    % hand back a slope of zero that looks like a measurement.
    % Scale the tolerance by the WAVENUMBER AXIS, not by kbar: in the degenerate case
    % kbar is itself ~0, so a self-relative tolerance can never trip.
    if (max(kbar) - min(kbar)) <= 1e-9 * max(kk)
        d.kbar = kbar;  d.fRidge = f(ok);  d.wRidge = (w(ok)/max(w(ok))).';
        return;
    end

    A    = [ones(nnz(ok),1), kbar(:)];

    wt_  = w(ok).' / max(w(ok));               % energy weights
    Aw   = A .* sqrt(wt_);
    yw   = om(:) .* sqrt(wt_);
    b    = Aw \ yw;
    r    = yw - Aw*b;
    omBar = sum(wt_ .* om(:)) / sum(wt_);

    d.slope   = b(2);
    d.slopeR2 = 1 - sum(r.^2) / max(sum(wt_ .* (om(:) - omBar).^2), eps);
    d.kbar    = kbar;
    d.fRidge  = f(ok);
    d.wRidge  = wt_.';

    cum = cumsum(wt_) / sum(wt_);
    lo  = find(cum >= 0.05, 1);  hi = find(cum >= 0.95, 1);
    d.fSupport = [d.fRidge(lo), d.fRidge(hi)];
end

% Author: Diellor Basha, 2026
