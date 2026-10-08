function ax = lambdaaxis(ax, Lambda, nTick)
% SHOW.LAMBDAAXIS  Label a mode-index y-axis with the spatial frequencies it corresponds to.
%
%   ord = rheome.show.lambdaaxis(ax, Lambda)          % after imagesc(ax, x, 1:K, Z(ord,:))
%
% WHY THIS EXISTS. An eigenbasis assembled per hemisphere is a CONCATENATION: Lambda runs
% 0 -> lambda_max for the left hemisphere and then restarts for the right. Passing that vector
% to imagesc as a y coordinate is wrong twice over -- imagesc uses only its first and last
% entries and spaces rows uniformly, so rows are drawn at positions they do not occupy, and the
% restart appears as a spurious horizontal discontinuity partway up the axis.
%
% Even after sorting, sqrt(lambda) is not uniformly spaced (Weyl: lambda_k ~ k, so sqrt(lambda_k)
% ~ sqrt(k)), so a uniform y coordinate still misplaces rows. The honest rendering is therefore
% ROW INDEX on the axis with sqrt(lambda) on the tick labels: no interpolation is implied and
% every row sits where it belongs.
%
% USAGE: sort your rows first --  [~,ord] = sort(Lambda);  imagesc(ax, x, 1:K, Z(ord,:));
% then call this to place the ticks.
%
% INPUTS:  ax handle;  Lambda [K x 1] eigenvalues (unsorted is fine);  nTick (default 8)
% OUTPUT:  ax
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(nTick), nTick = 8; end
    ls = sort(sqrt(double(Lambda(:))));
    K  = numel(ls);
    idx = unique(round(linspace(1, K, nTick)));
    set(ax, 'YTick', idx, 'YTickLabel', arrayfun(@(i) sprintf('%.0f', ls(i)), idx, 'UniformOutput', false));
end

% Author: Diellor Basha, 2026
