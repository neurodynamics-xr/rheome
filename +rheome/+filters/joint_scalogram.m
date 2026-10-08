function scal = joint_scalogram(C, Lambda, frame, f)
% FILTERS.JOINT_SCALOGRAM  Energy over SPATIAL SCALE x TEMPORAL FREQUENCY.
%
%   scal = rheome.filters.joint_scalogram(C, Lambda, frame, f)
%
% For each frame member m and each frequency bin,
%
%   E(m,omega) = || g_m(Lambda) .* C(:,omega) ||^2
%
% computed entirely in coefficients. Phi is M-orthonormal, so by Parseval this IS the spatial
% energy of the filtered field -- no synthesis to vertices, no [N x M x nOmega] volume.
%
% WHY THIS IS THE USEFUL ONE. A travelling structure satisfies omega = c*sqrt(lambda), so a
% constant-speed feature appears as a DIAGONAL RIDGE in the (sqrt(lambda), omega) plane. Reading
% its slope MEASURES the speed, which turns a speed-selective filter from an assumed velocity
% into one fitted to the data. The time-domain scalogram (rheome.filters.frame_scalogram) cannot show
% this: it has already integrated over frequency.
%
% INPUTS:
%   C       [K x nOmega] joint spectrum (from rheome.flow.joint)
%   Lambda  [K x 1] eigenvalues of the basis C is expressed in
%   frame   from rheome.filters.frame
%   f       [1 x nOmega] bin frequencies (Hz), for labelling
%
% OUTPUT (struct scal):
%   .energy  [M x nOmega]   per-member energy per frequency bin
%   .centers [1 x M]        member characteristic wavenumber sqrt(lambda)
%   .sigma   [1 x M]        member spatial scale (m), NaN for the scaling function
%   .f       [1 x nOmega]   frequencies
%   .slope   least-squares speed (m/s) from the ridge in (lambda,omega), and .slopeR2
%   .kbar    [1 x nOk] energy-weighted mean wavenumber per frequency (the ridge itself)
%   .fRidge  [1 x nOk] the frequencies it was evaluated at
%   .wRidge  [1 x nOk] normalised energy weight used in the slope fit
%   .fSupport [f5 f95] the frequency interval carrying 90%% of the energy
%
% ⚠ THE BAND MUST SPAN THE DIAGONAL. The slope is only recoverable if the retained frequency axis
% covers c*sqrt(lambda)/2pi across the wavenumbers present. Measured on synthetic diagonals with
% the energy-weighted fit (flow_joint_validate section 5):
%
%   c = 0.15 m/s -> +30.5%%   (diagonal clipped at the low-frequency edge of the grid)
%   c = 0.35 m/s ->  +6.1%%
%   c = 0.80 m/s ->  +2.2%%
%   c = 2.00 m/s ->  +3.3%%
%
% A few percent when the diagonal fits, biased UPWARD when it is clipped. A narrow band therefore
% cannot measure speed this way: an 8-13 Hz retention sees a 5 Hz slice of a diagonal that spans
% tens of Hz across cortical wavenumbers. Use a broadband retention for the speed estimate, and
% check .fSupport against the diagonal you expect before believing .slope.
%
% See also: rheome.filters.frame, rheome.filters.frame_scalogram, rheome.filters.frame_gains
%
% Author: Diellor Basha, 2026

    H = rheome.filters.frame_gains(frame, Lambda);        % [K x M]
    M = size(H,2);
    P = abs(C).^2;                                  % [K x nOmega]
    scal.energy = (H.^2).' * P;                     % [M x nOmega]  <- one GEMM
    scal.centers = frame.Centers;
    scal.f = double(f(:)).';

    % Member spatial scale. Prefer the frame's EXACT sigma = sqrt(2t); fall back to the centroid
    % only for a frame built before .Sigma existed. The centroid is a display summary contaminated
    % by where the eigenvalue axis was truncated, not the scale the member is matched to.
    if isfield(frame,'Sigma') && ~isempty(frame.Sigma)
        scal.sigma = frame.Sigma;
    else
        scal.sigma = sqrt(2) ./ max(frame.Centers, eps);
        if strcmpi(frame.Family,'mexhat'), scal.sigma(1) = NaN; end   % member 1 is the scaling fn
    end

    % ---- speed from the ridge, measured on the FULL (lambda,omega) energy ----
    % The frame-binned energy above is for display: with a handful of overlapping members the
    % wavenumber axis is far too coarse to fit a slope, and doing so returns nonsense. The diagonal
    % lives in (lambda,omega), so the ridge is taken there -- one energy-weighted mean wavenumber
    % per frequency bin -- and the frame is not involved.
    kk = sqrt(double(Lambda(:)));                   % [K x 1] wavenumbers
    w  = sum(P,1);                                  % [1 x nOmega] total energy per frequency
    ok = w > 0;
    if nnz(ok) > 2
        kbar = (kk.' * P(:,ok)) ./ w(ok);           % [1 x nOk] mean sqrt(lambda) per frequency
        om   = 2*pi*scal.f(ok);
        A    = [ones(nnz(ok),1), kbar(:)];

        % WEIGHTED by the energy at each frequency. Unweighted, every bin counts equally, so bins
        % carrying no signal -- outside the diagonal's span, where kbar is just the centroid of the
        % lambda axis -- dominate by sheer number and bias the slope badly (measured: +228% when the
        % diagonal filled only a third of the frequency axis). The ridge at a frequency is only as
        % informative as the energy there, so that is the weight.
        wt = w(ok).' / max(w(ok));
        Aw = A .* sqrt(wt);  yw = om(:) .* sqrt(wt);
        b  = Aw \ yw;
        r  = yw - Aw*b;
        omBar = sum(wt.*om(:))/sum(wt);
        scal.slope   = b(2);
        scal.slopeR2 = 1 - sum(r.^2)/max(sum(wt.*(om(:)-omBar).^2), eps);
        scal.kbar    = kbar;
        scal.fRidge  = scal.f(ok);
        scal.wRidge  = wt.';
        % effective support: where the energy actually is
        cum = cumsum(wt)/sum(wt);
        lo  = find(cum >= 0.05, 1);  hi = find(cum >= 0.95, 1);
        scal.fSupport = [scal.fRidge(lo), scal.fRidge(hi)];
    else
        scal.slope = NaN;  scal.slopeR2 = NaN;  scal.kbar = [];  scal.fRidge = [];
        scal.wRidge = [];  scal.fSupport = [NaN NaN];
    end
end

% Author: Diellor Basha, 2026
