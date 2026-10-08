% WHITEN  Noise-covariance regularization + whitener, shared by every +inverse solver.
%
%   iW = whiten(Cov, Method, NoiseReg, FourthMoment, nSamples)
%
% ⭐ SHARED ON PURPOSE. rheome.inverse.dirac and rheome.inverse.mne differ ONLY in the source basis, so any
% difference between their results must come from that and nothing else. Two copies of this
% function would make the comparison meaningless the first time one of them was edited.
%
% Port of bst_inverse_linear_2018 truncate_and_regularize_covariance (J.C. Mosher).
%
% Author: Diellor Basha, 2026

function iW = whiten(Cov, Method, NoiseReg, FourthMoment, nSamples)
% Port of bst_inverse_linear_2018 truncate_and_regularize_covariance (J.C. Mosher).
    Cov = (Cov + Cov')/2;
    [Un, Sn2] = svd(Cov, 'econ');
    Sn = sqrt(diag(Sn2));
    tol = length(Sn) * eps(single(Sn(1)));
    rk  = sum(Sn > tol);
    Un = Un(:, 1:rk);  Sn = Sn(1:rk);
    switch lower(Method)
        case 'none'
            iW = Un * diag(1./Sn) * Un';
        case 'median'
            Sn = max(Sn, median(Sn));
            iW = Un * diag(1./Sn) * Un';
        case 'diag'
            iW = diag(1 ./ sqrt(diag(Cov)));
        case 'reg'
            Ridge = mean(diag(Sn2)) * NoiseReg;               % Hamalainen mean-eigenvalue ridge
            iW = Un * diag(1 ./ sqrt(Sn.^2 + Ridge)) * Un';
        case 'shrink'
            [Cs, ~] = cov1para(Cov, FourthMoment, nSamples);
            [Un, Sn2] = svd(Cs, 'econ');  Sn = sqrt(diag(Sn2));
            tol = length(Sn) * eps(single(Sn(1)));  rk = sum(Sn > tol);
            Un = Un(:,1:rk);  Sn = Sn(1:rk);
            iW = Un * diag(1./Sn) * Un';
        otherwise
            error('inverse:whiten:method', 'Unknown NoiseMethod "%s".', Method);
    end
end

function [sNoiseCov, shrinkage] = cov1para(NoiseCov, FourthMoment, nSamples)
% Ledoit-Wolf shrinkage (verbatim from bst_inverse_linear_2018).
    n = size(NoiseCov,1);
    meanvar = mean(diag(NoiseCov));
    prior = meanvar * eye(n);
    phi = sum(sum(FourthMoment - NoiseCov.^2));
    gamma = norm(NoiseCov - prior, 'fro')^2;
    shrinkage = max(0, min(1, (phi/gamma)/nSamples));
    sNoiseCov = shrinkage*prior + (1-shrinkage)*NoiseCov;
end

function o = i_defaults(o, Def)
    fn = fieldnames(Def);
    for i = 1:numel(fn)
        if ~isfield(o, fn{i}) || isempty(o.(fn{i})), o.(fn{i}) = Def.(fn{i}); end
    end
end

% Author: Diellor Basha, 2026
