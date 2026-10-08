function q = qfactor(obj)
% QFACTOR  Centre wavenumber divided by half-power bandwidth.
%
%   q = qfactor(gfb)      [1 x M]
%
% Dimensionless, so this one transfers from cwtfilterbank unchanged -- a ratio is
% valid on any axis.
%
% See also: powerbw, centerWavenumbers
%
% Author: Diellor Basha, 2026
    bw = powerbw(obj);
    q  = centerWavenumbers(obj) ./ max((bw(:,2) - bw(:,1)).', eps);
end

% Author: Diellor Basha, 2026
