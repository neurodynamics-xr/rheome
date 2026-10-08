function [own, dominant] = crosstalk(vertexOperator, Gain, Membership)
% FLOW.CROSSTALK  What fraction of an ROI's measured flow actually originates in it.
%
%   [own, dominant] = rheome.flow.crosstalk(kc.vertexOperator, Gain, atlas.Membership)
%
%   own       [nROI x 1] fraction of the ROI's measured power that comes from inside it
%   dominant  [nROI x 1] index of the ROI the measured signal ACTUALLY comes from
%
% ⭐ THE ONLY HONEST BASIS FOR ROI ATTRIBUTION with a distributed inverse. R = K*Gain is the
% resolution matrix: row v says how a unit source at each location contributes to the measured
% flow at v. Aggregated over an ROI it says where that ROI's reading comes from. If most of it
% originates elsewhere, the attribution is leakage and NO NORMALISER REPAIRS IT.
%
% ⚠ DO NOT "CORRECT" BY DIVIDING BY SENSITIVITY. That was tried and it is wrong. Dividing
% measured energy by rheome.flow.sensitivity assumes a low-sensitivity vertex holds attenuated TRUE
% signal; it does not -- it holds leakage from elsewhere, and dividing amplifies exactly that.
% MEASURED on a reference subject: entorhinal R is 1.3% self-originating and parahippocampal R is
% 0.4%, yet dividing by sensitivity promoted both to the top of the alpha ranking, ahead of
% lateraloccipital R which is 67% self-originating. Use this as a MASK, not a divisor.
%
% ⚠ MOST ROIs DO NOT SURVIVE, and that is the finding rather than a failure. On the
% Desikan-Killiany parcellation with a 270-channel MEG minimum-norm inverse the MEDIAN own-ROI
% fraction is 6%; only 17 of 68 scouts reach 25%, together covering half the cortical area.
% Reporting a 68-row ROI table without saying which rows are resolvable overstates the
% instrument by a factor of four.
%
% ⚠ COST is one [nv x 3V] product per ROI -- about 5 s for 68 scouts at 20484 vertices and
% 270 channels. Computed once per anatomy and inverse, never per analysis.
%
% See also: rheome.flow.sensitivity, rheome.flow.curl, rheome.io.read.atlas
%
% Author: Diellor Basha, 2026

    K = double(vertexOperator);
    if size(K,2) ~= size(Gain,1)
        error('flow:crosstalk:size', ...
            'kernel has %d channels but the leadfield has %d rows.', size(K,2), size(Gain,1));
    end
    P  = double(Membership);
    nR = size(P, 1);
    own = zeros(nR,1);  dominant = zeros(nR,1);

    for j = 1:nR
        v = find(P(j,:));
        if isempty(v), own(j) = NaN; dominant(j) = NaN; continue; end
        p  = sum((K(v,:) * Gain).^2, 1);            % power per source COMPONENT
        pv = p(1:3:end) + p(2:3:end) + p(3:3:end);  % ...summed to a per-vertex source power
        tot = sum(pv);
        if ~(tot > 0), own(j) = NaN; dominant(j) = NaN; continue; end
        pv = pv(:) / tot;
        own(j) = sum(pv(v));
        [~, dominant(j)] = max(P * pv);
    end
end

% Author: Diellor Basha, 2026
