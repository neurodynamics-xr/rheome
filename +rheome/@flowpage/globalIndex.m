function gi = globalIndex(obj, lambdaCut)
% GLOBALINDEX  Fraction of spectral energy below a lambda cut, per frame.
%
%   gi = globalIndex(fp)              % cut at the page's ENERGY-median lambda
%   gi = globalIndex(fp, lambdaCut)   % fixed cut, comparable across pages
%
% The running global/local readout: 1 means all the energy sits in globally distributed low
% spatial modes, 0 means it is all in spatially confined high modes. Collapses the spectrum
% to one number so it can be watched evolving on a timeline; scaleEnergy keeps the shape.
%
% ⚠ LAMBDA IS NOT SORTED. rheome.flow.context assembles the basis BLOCK-DIAGONALLY per hemisphere,
% so mode index is not lambda order and a cumulative sum over modes is not a cumulative sum
% over lambda. Everything spectral here sorts first. (A weighted mean such as the centroid is
% order-independent and unaffected; a cumulative or a line plot is not.)
%
% ⚠ THE DEFAULT CUT IS THE PAGE'S OWN ENERGY MEDIAN, not the median of the lambda axis.
% MEASURED: with the axis median the index sat at 0.94-0.96 across 4000 frames of resting
% alpha -- saturated, and useless for watching anything change, because alpha vorticity is
% overwhelmingly coarser than the axis midpoint (~40 mm). Centring the cut on the page puts
% the index near 0.5 so its VARIATION is the signal. That makes it page-relative: pass an
% explicit lambdaCut when comparing pages or subjects.
%
% See also: scaleEnergy, centroid
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda(:);
    if nargin < 2 || isempty(lambdaCut)
        Pall = sum(abs(double(obj.Coefficients)).^2, 2);      % energy per mode over the page
        [lamS, ord] = sort(lam);                             % LAMBDA order, not mode order
        cum  = cumsum(Pall(ord)) / max(sum(Pall), realmin);
        iCut = find(cum >= 0.5, 1);
        if isempty(iCut), iCut = numel(lamS); end
        lambdaCut = lamS(iCut);
    end
    low = lam <= lambdaCut;

    P  = abs(double(obj.Coefficients)).^2;        % [Ks x nT]
    tot = sum(P, 1);
    gi  = sum(P(low, :), 1) ./ max(tot, realmin);
    gi(tot <= 0) = NaN;                            % no energy -> undefined, not 0
end

% Author: Diellor Basha, 2026
