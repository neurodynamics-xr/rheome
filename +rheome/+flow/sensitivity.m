function s = sensitivity(vertexOperator, noiseCov)
% FLOW.SENSITIVITY  What a flow kernel produces from NOISE, per vertex.
%
%   s = rheome.flow.sensitivity(kc.vertexOperator, noiseCov)     -> [V x 1]
%
%   s_v = k_v * C_noise * k_v'      (the variance the kernel yields at v from that noise)
%
% ⭐ THE NORMALISER THAT MAKES REGIONS COMPARABLE. A measured flow energy is not interpretable
% on its own, because the kernel's sensitivity varies enormously across the surface -- deep and
% medial vertices are poorly constrained and a handful of vertices dominate. Dividing measured
% energy by this makes the result a RATIO to what noise alone would produce there: dimensionless,
% comparable between regions of any size or depth, and directly readable as an SNR.
%
% ⚠ THIS IS A BUG FIX, NOT A REFINEMENT. Measured on the AnphySleep EEG kernel: the top 1% of
% vertices carry 15% of the total squared row norm, and the smallest scouts (entorhinal,
% frontalpole, transversetemporal) are exactly the ones with the largest kernel norm PER UNIT
% AREA. Reporting area-normalised ROI energy without this returned 'transversetemporal R' as the
% peak region for EVERY band in EVERY sleep stage -- a property of the inverse, not the data.
% The same mechanism inflated the deep and medial regions in the resting MEG phase-locking map.
%
% ⚠ REDUCE THE NORMALISER THE SAME WAY AS THE MEASUREMENT. To normalise an ROI mean, reduce s
% over the ROI with the SAME area weighting, then divide:
%     roiS = (Membership * (w .* s)) ./ areas;      roiRatio = roiE ./ roiS;
% Normalising a per-area mean by a per-vertex sensitivity does not cancel the bias.
%
% ⚠ WITH AN IDENTITY COVARIANCE this is just the squared row norm -- the kernel's gain to
% isotropic input. That is a reasonable fallback when no noise covariance is available, but the
% measured covariance is the better null because it carries the real spatial structure of the
% noise.
%
% INPUTS:
%   vertexOperator  [V x C] sensors -> per-vertex flow quantity (kc.vertexOperator)
%   noiseCov        [C x C] sensor noise covariance, or [] for identity
%
% See also: rheome.flow.curl, rheome.flowfeatures, rheome.io.read.noisecov
%
% Author: Diellor Basha, 2026

    K = double(vertexOperator);
    C = size(K, 2);
    if nargin < 2 || isempty(noiseCov)
        s = sum(K.^2, 2);
        return;
    end
    Cn = double(noiseCov);
    if ~isequal(size(Cn), [C C])
        error('flow:sensitivity:size', ...
            'noiseCov is %s but the kernel has %d channels.', mat2str(size(Cn)), C);
    end
    Cn = (Cn + Cn.') / 2;                       % symmetrise; a covariance read from disk may not be
    s = sum((K * Cn) .* K, 2);
    s = max(s, 0);                              % a PSD-to-rounding covariance can give -1e-30
end

% Author: Diellor Basha, 2026
