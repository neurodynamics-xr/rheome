function varargout = geodesic(varargin)
% GEOM.GEODESIC  Deprecated: renamed rheome.geom.measure.heatDistance (the heat-method geodesic distance).
%
%   [distance, solver] = rheome.geom.geodesic(S, sourceVertices, solver)
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to rheome.geom.measure.heatDistance
% with the same inputs and outputs, and warns once per session (rheome:deprecated). Call rheome.geom.measure.heatDistance.
%
% See also: rheome.geom.measure.heatDistance
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.geodesic is deprecated: call rheome.geom.measure.heatDistance (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.measure.heatDistance(varargin{:});
end
