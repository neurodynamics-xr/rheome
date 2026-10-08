function [ok, why] = iscompatible(obj, ax, tol)
% ISCOMPATIBLE  Do this bank and another representation live on the same axes?
%
%   [ok, why] = iscompatible(jfb, ax)
%
% ⚠ .half and .boundary are checked because a mismatch there is silent and fatal. A
% positive-half (analytic, complex) representation and a full-axis (real) one cannot be
% combined, nor can a ring-boundary and a path-boundary one, whose implied time
% extensions differ. And two runs over DIFFERENT bands can retain the SAME number of
% bins, so a size check passes them both.
%
% ⭐ THE COUNTS AND CONVENTIONS ALSO LIVE IN A TYPE NOW: jointdomain(jfb) returns the joint plane as a
% rheome.domain.product and rheome.domain.same compares its factors and its conventions, so that half of this check
% exists once rather than per class. ⚠ The VALUE comparison below cannot move there: a domain carries
% no coordinates, so no id comparison can catch the same mode count over different anatomy. Both are
% needed and neither subsumes the other.
%
% See also: axes, jointdomain, rheome.domain.same
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(tol), tol = 1e-10; end
    mine = axes(obj);
    why  = {};

    if isfield(ax,'lambda') && ~isempty(ax.lambda)
        if numel(ax.lambda) ~= numel(mine.lambda)
            why{end+1} = sprintf('lambda: %d vs %d modes', ...
                numel(mine.lambda), numel(ax.lambda));
        elseif max(abs(mine.lambda - ax.lambda(:))) > tol*max(max(mine.lambda), eps)
            why{end+1} = 'lambda: same count, different values (different anatomy or basis)';
        end
    end
    if isfield(ax,'f') && ~isempty(ax.f)
        if numel(ax.f) ~= numel(mine.f)
            why{end+1} = sprintf('frequency: %d vs %d bins -- different bands retained?', ...
                numel(mine.f), numel(ax.f));
        elseif max(abs(mine.f - ax.f(:).')) > tol*max(max(mine.f), eps)
            why{end+1} = 'frequency: same count, different values';
        end
    end
    if isfield(ax,'half') && ~strcmpi(mine.half, ax.half)
        why{end+1} = sprintf('half: ''%s'' vs ''%s'' -- one is analytic/complex, the other real', ...
            mine.half, ax.half);
    end
    if isfield(ax,'boundary') && ~strcmpi(mine.boundary, ax.boundary)
        why{end+1} = sprintf('boundary: ''%s'' vs ''%s'' -- different implied time extension', ...
            mine.boundary, ax.boundary);
    end
    ok = isempty(why);
end

% Author: Diellor Basha, 2026
