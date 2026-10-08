function fit = pde_fit(coeff, opts)
% DYNAMICS.PDE_FIT  Fit the per-mode temporal dynamics to candidate PDEs (the model-based law).
%
%   fit = rheome.dynamics.pde_fit(coeff [,opts])
%
% For each spatial mode k (eigenvalue lambda_k), fits its coefficient trajectory c_k(t) to two
% linear generators and extracts the effective physical constants:
%   diffusion    c_k(t) = c_k(0) exp(-lambda_k D t)       -> effective diffusivity D
%   damped wave  c_k(t) = exp(-gamma t) cos(sqrt(lambda_k) c t + phi)  -> wave speed c, damping gamma
% Each mode is classified (oscillatory -> wave; monotone decay -> diffusion) and the constants are
% aggregated (robust median) across modes. This is the model-based counterpart to the data-driven
% DMD/dispersion; the two should agree on the wave speed.
%
% INPUTS:
%   coeff  from rheome.dynamics.mode_coefficients (.C [K x T], .dt, .lambda [K x 1])
%   opts.minCycles  minimum oscillation cycles to accept a wave mode (default 1.5)
% OUTPUT (struct fit):
%   .diffusivity  robust D over diffusion modes           (position^2/s)
%   .waveSpeed    robust c over wave modes                 (position/s)
%   .damping      robust gamma over wave modes             (1/s)
%   .modelType    'wave' | 'diffusion' (majority)
%   .perMode      struct array: .lambda .frequency .speed .damping .diffusivity .isWave
%
% See also: rheome.dynamics.dmd, rheome.dynamics.dispersion
%
% Author: Diellor Basha, 2026

    if nargin < 2, opts = struct(); end
    if ~isfield(opts,'minCycles') || isempty(opts.minCycles), opts.minCycles = 1.5; end
    C = coeff.C;  dt = coeff.dt;  lambda = coeff.lambda(:);
    [K, T] = size(C);  duration = (T-1)*dt;
    freqAxis = (0:T-1)/(T*dt);

    speedList = []; dampingList = []; diffusivityList = []; isWaveCount = 0;
    perMode = struct('lambda',{},'frequency',{},'speed',{},'damping',{},'diffusivity',{},'isWave',{});
    for k = 1:K
        x = C(k,:) - mean(C(k,:));
        if std(x) < 1e-12 || lambda(k) <= 0, continue; end
        % dominant temporal frequency
        spec = abs(fft(x)).^2;  spec(1) = 0;  half = 1:floor(T/2);
        [~, iPk] = max(spec(half));  peakFreq = freqAxis(iPk);
        numCycles = peakFreq * duration;
        envelope = abs(hilbert(x));                                   % amplitude envelope
        trim = round(0.05*T);  inner = (trim+1):(T-trim);            % drop Hilbert edge artifacts
        good = inner(envelope(inner) > 0.15*max(envelope(inner)));    % high-SNR portion only
        dampingSlope = polyfit((good-1)*dt, log(envelope(good)), 1);  % log-env vs t
        gamma = -dampingSlope(1);
        isWave = numCycles >= opts.minCycles;
        pm = struct('lambda',lambda(k),'frequency',peakFreq,'speed',NaN,'damping',gamma, ...
                    'diffusivity',NaN,'isWave',isWave);
        if isWave
            pm.speed = 2*pi*peakFreq / sqrt(lambda(k));               % omega = c*sqrt(lambda)
            speedList(end+1) = pm.speed;   dampingList(end+1) = gamma;  %#ok<AGROW>
            isWaveCount = isWaveCount + 1;
        else
            % monotone-ish decay: diffusivity from log|c_k| slope,  slope = -lambda_k * D
            mag = abs(C(k,:));  m2 = mag > 0.05*max(mag);
            decaySlope = polyfit((find(m2)-1)*dt, log(mag(m2)), 1);
            pm.diffusivity = -decaySlope(1) / lambda(k);
            diffusivityList(end+1) = pm.diffusivity;                  %#ok<AGROW>
        end
        perMode(end+1) = pm;                                          %#ok<AGROW>
    end

    fit = struct();
    fit.waveSpeed    = median(speedList(isfinite(speedList)));
    fit.damping      = median(dampingList(isfinite(dampingList)));
    fit.diffusivity  = median(diffusivityList(isfinite(diffusivityList)));
    fit.modelType    = i_pick(isWaveCount >= numel(perMode)-isWaveCount, 'wave', 'diffusion');
    fit.perMode      = perMode;
end

function s = i_pick(cond, a, b)
    if cond, s = a; else, s = b; end
end

% Author: Diellor Basha, 2026
