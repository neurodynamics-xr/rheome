function relation = dispersion(dmdResult, lambdaPerMode)
% DYNAMICS.DISPERSION  Empirical dispersion relation omega(lambda) from a DMD result.
%
%   relation = rheome.dynamics.dispersion(dmdResult, lambdaPerMode)
%
% For each DMD mode, pairs its temporal angular frequency omega = 2*pi*|f| with its dominant SPATIAL
% scale (the eigenvalue lambda at which the mode's spatial power concentrates), giving the empirical
% dispersion relation of the cortical dynamics. A linear ridge omega = c*sqrt(lambda) is a traveling
% WAVE of phase speed c; omega proportional to lambda is DIFFUSION. The model with the better fit is
% reported, with the wave-speed estimate.
%
% INPUTS:
%   dmdResult      output of rheome.dynamics.dmd (.modes [K x r], .frequency [r x 1])
%   lambdaPerMode  [K x 1] spatial eigenvalue of each basis mode (coeff.lambda)
% OUTPUT (struct relation):
%   .omega .lambda   the (omega, lambda) points used
%   .speedFit        fitted wave speed c from omega = c*sqrt(lambda)   (position-units / s)
%   .modelType       'wave' | 'diffusion' (whichever fits better)
%   .rWave .rDiff    fit quality (R^2) of each model
%
% See also: rheome.dynamics.dmd, rheome.filters.travwave
%
% Author: Diellor Basha, 2026

    modePower = abs(dmdResult.modes).^2;                                  % [K x r]
    dominantLambda = (lambdaPerMode(:)' * modePower)' ./ max(sum(modePower, 1)', eps);  % power-weighted
    omega = 2*pi*abs(dmdResult.frequency);
    keep = omega > 0 & isfinite(dominantLambda) & dominantLambda > 0;
    omega = omega(keep);  dominantLambda = dominantLambda(keep);

    sqrtLambda = sqrt(dominantLambda);
    % wave: omega = c*sqrt(lambda) (through origin);  diffusion: omega = D*lambda (through origin)
    speedFit = (sqrtLambda' * omega) / max(sqrtLambda' * sqrtLambda, eps);
    diffusFit = (dominantLambda' * omega) / max(dominantLambda' * dominantLambda, eps);
    totalVar = max(sum((omega - mean(omega)).^2), eps);
    rWave = 1 - sum((omega - speedFit*sqrtLambda).^2) / totalVar;
    rDiff = 1 - sum((omega - diffusFit*dominantLambda).^2) / totalVar;

    relation = struct('omega', omega, 'lambda', dominantLambda, ...
        'speedFit', speedFit, 'diffusivityFit', diffusFit, ...
        'modelType', i_pick(rWave >= rDiff, 'wave', 'diffusion'), ...
        'rWave', rWave, 'rDiff', rDiff);
end

function s = i_pick(cond, a, b)
    if cond, s = a; else, s = b; end
end

% Author: Diellor Basha, 2026
