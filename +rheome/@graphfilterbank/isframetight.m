function tf = isframetight(obj, tol)
% ISFRAMETIGHT  Does this bank reconstruct without a dual?
%   tf = isframetight(gfb)        tol = 1e-6 on B/A - 1
%
% See also: framebounds
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(tol), tol = 1e-6; end
    tf = abs(framebounds(obj).Tightness - 1) <= tol;
end

% Author: Diellor Basha, 2026
