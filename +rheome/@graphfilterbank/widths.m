function s = widths(obj)
% WIDTHS  Member spatial width sigma = sqrt(2t). EXACT, not fitted: a mexhat member
% peaks at lambda = 1/t and a Gaussian of width sigma is exp(-lambda*sigma^2/2), so
% t = sigma^2/2 identically.
%
%   s = widths(gfb)      [1 x M], in the operator's length units
%
% ⭐ REPORT THIS, not centerWavenumbers, for member scale. The centroid is a display
% summary contaminated by where the spectrum was truncated.
%
% See also: scales, centerWavenumbers
%
% Author: Diellor Basha, 2026
    s = sqrt(2 * obj.T_);
end

% Author: Diellor Basha, 2026
