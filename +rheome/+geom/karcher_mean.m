function [centerVertex, centerPos, dispersion] = karcher_mean(S, points, opts)
% GEOM.KARCHER_MEAN  Intrinsic (geodesic Fréchet) mean of points on the surface.
%
%   [centerVertex, centerPos, dispersion] = rheome.geom.karcher_mean(S, points)
%   [...] = rheome.geom.karcher_mean(S, points, opts)
%
% The surface point minimising sum_i w_i * d(p, x_i)^power (geodesic distance) -- the discrete
% Fréchet mean. power=2 is the Karcher mean; power=1 the geometric median. On a mesh this is an
% argmin over vertices of a per-vertex objective built from one rheome.geom.geodesic field per source
% point (the cached solver makes the repeated queries cheap). Faithful to nxr-compute's
% definition (compute.h: "minimizes sum d(p, source_i)^p"); no log-map / transport port needed.
% Unlike the ambient (Euclidean) mean, the result stays ON the sheet.
%
% INPUTS:
%   S      surface struct (.Vertices [nV x 3], .Faces)
%   points [n x 3] positions (snapped to nearest vertex) OR [n x 1] vertex indices
%   opts   .weights [n x 1] (default ones)  .power (default 2)  .solver (reuse rheome.geom.geodesic)
%
% OUTPUT:
%   centerVertex  index of the minimiser vertex
%   centerPos     [1 x 3] its position
%   dispersion    ( weighted-mean geodesic d^power )^(1/power) = intrinsic spread (metres)
%
% See also: rheome.geom.geodesic
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct; end
    if ~isfield(opts,'power') || isempty(opts.power), opts.power = 2; end
    if size(points,2) == 3
        sourceVertices = zeros(size(points,1),1);
        for i = 1:size(points,1)
            [~, sourceVertices(i)] = min(sum((S.Vertices - points(i,:)).^2, 2));
        end
    else
        sourceVertices = points(:);
    end
    n = numel(sourceVertices);
    if ~isfield(opts,'weights') || isempty(opts.weights), opts.weights = ones(n,1); end
    weights = opts.weights(:);

    if isfield(opts,'solver'), solver = opts.solver; else, solver = []; end
    objective = zeros(size(S.Vertices,1), 1);
    for i = 1:n
        [distance, solver] = rheome.geom.geodesic(S, sourceVertices(i), solver);
        objective = objective + weights(i) * distance.^opts.power;
    end
    [minObjective, centerVertex] = min(objective);
    centerPos  = S.Vertices(centerVertex, :);
    dispersion = (minObjective / sum(weights))^(1/opts.power);
end

% Author: Diellor Basha, 2026
