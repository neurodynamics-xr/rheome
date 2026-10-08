function tf = isframetight(obj, tol)
% ISFRAMETIGHT  True if B/A - 1 <= tol (default 1e-9).
%   tf = isframetight(tfb [, tol])
% Author: Diellor Basha, 2026
    if nargin < 2, tol = 1e-9; end
    b = framebounds(obj);
    tf = b.Tightness - 1 <= tol;
end
% Author: Diellor Basha, 2026
