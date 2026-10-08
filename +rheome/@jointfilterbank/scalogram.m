function E = scalogram(obj, C)
% SCALOGRAM  Energy per member per temporal-frequency bin.
%
%   E = scalogram(jfb, C)      C [K x nOmega] -> E [Nf x nOmega]
%
%   E(i, omega) = || W_i(:,omega) .* C(:,omega) ||^2
%
% Computed entirely in joint coefficients and ONE MEMBER AT A TIME, so it works below the
% byte guard where wt does not: if you want the picture rather than the coefficients,
% this is the call.
%
% ⚠ Do NOT fit a dispersion slope to this. With a handful of overlapping members the
% wavenumber axis is far too coarse; use dispersion, which works on the full
% (lambda, omega) energy.
%
% See also: dispersion, wt, framebounds
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda;  om = obj.Omega;
    K = numel(lam);  nO = numel(om);
    if ~isequal(size(C), [K nO])
        error('jointfilterbank:size', ...
            'C is %s but the bank grid is %s.', mat2str(size(C)), mat2str([K nO]));
    end
    E = zeros(obj.NumMembers, nO);
    for m = 1:obj.NumMembers
        E(m, :) = sum(abs(C .* jfb_member(obj, m, lam, om)).^2, 1);
    end
end

% Author: Diellor Basha, 2026
