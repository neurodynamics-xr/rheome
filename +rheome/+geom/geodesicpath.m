function varargout = geodesicpath(varargin)
% GEOM.GEODESICPATH  Deprecated: renamed rheome.geom.draw.line (a geodesic between two vertices, on the faces).
%
%   [P, info] = rheome.geom.geodesicpath(S, from, to, solver)   % solver: ignored (see rheome.geom.draw.line), comes back as info.solver
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to rheome.geom.draw.line
% with the same inputs and outputs, and warns once per session (rheome:deprecated). Call rheome.geom.draw.line.
%
% See also: rheome.geom.draw.line
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.geodesicpath is deprecated: call rheome.geom.draw.line (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.draw.line(varargin{:});
end
