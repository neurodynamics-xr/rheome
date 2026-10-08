function result = dmd(C, opts)
% DYNAMICS.DMD  Reduced-rank Dynamic Mode Decomposition of a coefficient time series.
%
%   result = rheome.dynamics.dmd(C [,opts])
%
% Fits the linear propagator A in C(:,t+1) = A C(:,t) via a truncated-SVD (reduced-rank) DMD, and
% returns its eigenvalues/modes -- the effective dynamical operator of the (reduced-order) cortical
% activity. No transport/advection assumption: a traveling wave appears as a complex mode with a
% spatial phase gradient. Each DMD eigenvalue mu gives a temporal frequency and growth/decay rate.
%
% INPUTS:
%   C     [K x T] coefficient time series (real or complex; from rheome.dynamics.mode_coefficients)
%   opts.dt      sampling interval (s), default 1
%   opts.rank    SVD truncation rank; default = 99%-energy knee (>=2)
%   opts.delays  time-delay (Hankel) embedding depth, default 1 (none). Use >1 to resolve REAL
%                standing oscillations into their +/-omega pairs (needed for real-valued MEG).
% OUTPUT (struct result):
%   .eigenvalues [r x 1]  DMD eigenvalues mu_k
%   .frequency   [r x 1]  temporal frequency f_k = angle(mu_k)/(2*pi*dt)   (Hz)
%   .growth      [r x 1]  growth/decay g_k = log|mu_k|/dt                  (1/s; <0 = decaying)
%   .modes       [K x r]  DMD modes in coefficient space (reconstruct to cortex via the basis)
%   .amplitude   [r x 1]  mode amplitudes (projection of the initial state)
%   .rank        scalar   truncation rank used
%
% See also: rheome.dynamics.mode_coefficients, rheome.dynamics.dispersion, rheome.dynamics.pde_fit
%
% Author: Diellor Basha, 2026

    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'dt')     || isempty(opts.dt),     opts.dt = 1;     end
    if ~isfield(opts, 'delays') || isempty(opts.delays), opts.delays = 1; end

    numChannels = size(C, 1);
    % time-delay (Hankel) embedding: stack the signal with its shifts so real standing oscillations
    % are resolvable into +/-omega pairs.
    if opts.delays > 1
        d = opts.delays;  Tcols = size(C,2) - d + 1;
        Cembed = zeros(numChannels*d, Tcols);
        for s = 1:d
            Cembed((s-1)*numChannels + (1:numChannels), :) = C(:, s:s+Tcols-1);
        end
        C = Cembed;
    end

    snapPast   = C(:, 1:end-1);
    snapFuture = C(:, 2:end);
    [leftSing, singVals, rightSing] = svd(snapPast, 'econ');
    singular = diag(singVals);

    if ~isfield(opts, 'rank') || isempty(opts.rank)
        rank = find(cumsum(singular)/sum(singular) >= 0.99, 1);   % 99%-energy knee
        rank = max(2, min(rank, numel(singular)));
    else
        rank = min(opts.rank, numel(singular));
    end

    Ur = leftSing(:, 1:rank);  Sr = singVals(1:rank, 1:rank);  Vr = rightSing(:, 1:rank);
    reducedPropagator = Ur' * snapFuture * Vr / Sr;
    [eigVec, eigMat] = eig(reducedPropagator);
    mu = diag(eigMat);
    dmdModes = snapFuture * Vr / Sr * eigVec;                     % [K*delays x r] (embedded space)
    amplitude = dmdModes \ C(:, 1);                               % initial-state projection
    spatialModes = dmdModes(1:numChannels, :);                    % first delay block = spatial mode

    result = struct('eigenvalues', mu, ...
                    'frequency', angle(mu)/(2*pi*opts.dt), ...
                    'growth',    log(abs(mu))/opts.dt, ...
                    'modes',     spatialModes, ...
                    'amplitude', amplitude, ...
                    'rank',      rank);
end

% Author: Diellor Basha, 2026
