function [hper, hap] = gains(P, Pap)
% SPECTRAL.GAINS  Complementary periodic / aperiodic gains from a fitted background.
%
%   [hper, hap] = rheome.spectral.gains(P, Pap)
%
% Returns the pair of real gains that split a field into periodic and aperiodic parts:
%
%   hap  = sqrt( Pap ./ P )                 aperiodic share of the AMPLITUDE
%   hper = sqrt( max(0, 1 - Pap./P) )       periodic share  = the residual
%
% so that hper.^2 + hap.^2 == 1 pointwise, and applying them to a complex coefficient array
% partitions POWER exactly: |C_per|^2 + |C_ap|^2 == |C|^2.
%
% WHY THE RESIDUAL, NOT THE FITTED PEAKS. specparam is additive in LOG power, hence multiplicative
% in linear power, so sqrt(peaks./P) does not complement sqrt(Pap./P) and the two would not
% partition. Defining the periodic part as whatever exceeds the background makes the partition
% exact by construction, and makes the split independent of how well the peaks were fitted.
%
% CLIPPING. Where the fit overshoots (Pap > P) the residual would be negative; it is clipped to
% zero, which sets hper = 0 and hap = 1 there. This is renormalized so the identity still holds
% exactly. Clipped bins are reported.
%
% INPUTS:
%   P    [nF x nS] linear power (>= 0)
%   Pap  [nF x nS] fitted aperiodic power, same size (from rheome.spectral.aperiodic)
%
% OUTPUT:
%   hper, hap  [nF x nS] real, in [0,1], with hper.^2 + hap.^2 = 1
%
% See also: rheome.spectral.aperiodic, rheome.spectral.decompose
%
% Author: Diellor Basha, 2026

    if ~isequal(size(P), size(Pap))
        error('spectral:gains:size', 'P is %s but Pap is %s.', mat2str(size(P)), mat2str(size(Pap)));
    end
    P    = max(double(P), realmin);
    frac = min(max(double(Pap) ./ P, 0), 1);      % aperiodic share of POWER, clipped to [0,1]
    hap  = sqrt(frac);
    hper = sqrt(1 - frac);
end

% Author: Diellor Basha, 2026
