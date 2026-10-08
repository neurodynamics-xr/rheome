function [U, gt] = fbd_rotating(S, degrees, Omega, tvec)
% FBD_ROTATING  A rigidly rotating pattern, and the ridge it must produce.
%
%   [U, gt] = fbd_rotating(S, degrees, Omega, tvec)
%
% Superposes sectoral harmonics of the given degrees, all rotating rigidly at Omega:
%
%   u(v,t) = sum_l sin^l(theta) * cos( l*(phi - Omega*t) )
%
% ⭐ THE RIDGE IS ANALYTIC. Degree l contributes temporal frequency omega_l = l*Omega at
% wavenumber k_l = sqrt(l(l+1))/R, so the energy lies on a locus in (k, omega) whose
% least-squares slope is the speed dispersion() must recover. gt.slope is that slope,
% computed from the planted pairs -- never fitted to the output.
%
% OUTPUT: U [nV x nT];  gt.k, gt.omega [1 x nL], gt.slope, gt.degrees
%
% Author: Diellor Basha, 2026

    degrees = degrees(:).';
    nT = numel(tvec);
    U  = zeros(S.nV, nT);
    for l = degrees
        for it = 1:nT
            U(:, it) = U(:, it) + fbd_sectoral(S, l, -l*Omega*tvec(it));
        end
    end

    k     = sqrt(degrees .* (degrees + 1)) / S.R;
    omega = degrees * Omega;
    if numel(degrees) >= 2
        p = polyfit(k, omega, 1);         % the locus the planted content lies on
        slope = p(1);
    else
        % One degree is one point: a line through it is not determined, but the PHASE
        % SPEED omega/k is, and that is the physically meaningful answer.
        slope = omega / k;
    end
    gt = struct('k', k, 'omega', omega, 'slope', slope, 'degrees', degrees);
end

% Author: Diellor Basha, 2026
