function varargout = flipOut(varargin)
% GEOM.FLIPOUT  Deprecated: renamed rheome.geom.draw.flipOut (the FlipOut geodesic between two vertices).
%
%   [P, info] = rheome.geom.flipOut(S, from, to, mesh)          % mesh = rheome.geom.draw.flipMesh(S) (formerly info.cache)
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to rheome.geom.draw.flipOut
% with the same inputs and outputs, and warns once per session (rheome:deprecated). Call rheome.geom.draw.flipOut.
%
% See also: rheome.geom.draw.flipOut
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.flipOut is deprecated: call rheome.geom.draw.flipOut (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.draw.flipOut(varargin{:});
end
