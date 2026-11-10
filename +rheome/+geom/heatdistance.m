function varargout = heatdistance(varargin)
% GEOM.HEATDISTANCE  Deprecated: renamed rheome.geom.measure.heatDistanceEigen (the heat-method distance in the eigenbasis).
%
%   phi = rheome.geom.heatdistance(basis, fg, seeds, opts)
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to rheome.geom.measure.heatDistanceEigen
% with the same inputs and outputs, and warns once per session (rheome:deprecated). Call rheome.geom.measure.heatDistanceEigen.
%
% See also: rheome.geom.measure.heatDistanceEigen
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.heatdistance is deprecated: call rheome.geom.measure.heatDistanceEigen (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.measure.heatDistanceEigen(varargin{:});
end
