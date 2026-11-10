function varargout = tilemean(varargin)
% GEOM.TILEMEAN  Deprecated: renamed rheome.geom.measure.tileMean (the per-tile mean or median of a vertex map).
%
%   Y = rheome.geom.tilemean(G, X, a, stat)
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to rheome.geom.measure.tileMean
% with the same inputs and outputs, and warns once per session (rheome:deprecated). Call rheome.geom.measure.tileMean.
%
% See also: rheome.geom.measure.tileMean
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.tilemean is deprecated: call rheome.geom.measure.tileMean (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.measure.tileMean(varargin{:});
end
