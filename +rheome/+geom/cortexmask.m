function varargout = cortexmask(varargin)
% GEOM.CORTEXMASK  Deprecated: renamed rheome.geom.select.cortex (the vertices that are cortex).
%
%   c = rheome.geom.cortexmask(src, nV)
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to rheome.geom.select.cortex
% with the same inputs and outputs, and warns once per session (rheome:deprecated). Call rheome.geom.select.cortex.
%
% See also: rheome.geom.select.cortex
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.cortexmask is deprecated: call rheome.geom.select.cortex (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.select.cortex(varargin{:});
end
