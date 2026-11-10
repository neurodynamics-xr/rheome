function varargout = karcher_mean(varargin)
% GEOM.KARCHER_MEAN  Deprecated: renamed rheome.geom.measure.karcherMean (the intrinsic mean of points on the surface).
%
%   [centerVertex, centerPos, dispersion] = rheome.geom.karcher_mean(S, points, opts)
%
% Moved 2026-11-10 when rheome.geom's surface functions were grouped by job into measure (numbers), select (a set of
% vertices) and draw (a curve on the surface), the Editor's and nxr's names. This name forwards to
% rheome.geom.measure.karcherMean with the same inputs and outputs, and warns once per session (rheome:deprecated).
% Call rheome.geom.measure.karcherMean.
%
% See also: rheome.geom.measure.karcherMean
%
% Author: Diellor Basha, 2026

    persistent warned
    if isempty(warned)
        warned = true;
        warning('rheome:deprecated', 'rheome.geom.karcher_mean is deprecated: call rheome.geom.measure.karcherMean (same inputs and outputs).');
    end
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = rheome.geom.measure.karcherMean(varargin{:});
end
