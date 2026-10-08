function T = rows(analysis, metric, value, unit, band)
% SCALE.ROWS  Long-format metric rows: one row per (analysis, metric, band).
%
%   T = rheome.scale.rows("resolution", ["r50_median" "ple_median"], [52 23], "mm")
%   T = rheome.scale.rows(analysis, metric, value, unit, band)   % band defaults to ""
%
% The per-subject table every scale.measure_* returns; rheome.scale.run adds subject and cohort, and
% rheome.scale.reduce reads only this shape.
%
% See also: rheome.scale.run, rheome.scale.reduce
%
% Author: Diellor Basha, 2026

    if nargin < 5 || isempty(band), band = ""; end
    metric = string(metric(:));  n = numel(metric);
    value = double(value(:));  unit = i_n(string(unit), n);  band = i_n(string(band), n);
    T = table(repmat(string(analysis), n, 1), metric, band, value, unit, ...
              'VariableNames', {'analysis','metric','band','value','unit'});
end

function s = i_n(s, n)
    s = s(:);  if isscalar(s), s = repmat(s, n, 1); end
end

% Author: Diellor Basha, 2026
