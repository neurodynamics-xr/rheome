function theta = sel_threshold(thr, V)
% SEL_THRESHOLD  Resolve a threshold: a number, "pNN", or "a*median", over the scope's values.
% Author: Diellor Basha, 2026
    if isnumeric(thr), theta = double(thr); return; end
    t = char(string(thr));
    v = V(:);  v = v(isfinite(v));
    tok = regexp(t, '^p(\d+(\.\d+)?)$', 'tokens', 'once');
    if ~isempty(tok), theta = i_prctile(v, str2double(tok{1})); return; end
    tok = regexp(t, '^(\d+(\.\d+)?)\*median$', 'tokens', 'once');
    if ~isempty(tok), theta = str2double(tok{1}) * median(v); return; end
    error('select:threshold', 'Threshold must be a number, "pNN" or "a*median"; got %s.', t);
end

function y = i_prctile(x, p)
% MATLAB's prctile definition (plotting positions 100*(i-0.5)/n, linear interpolation,
% clamped at the ends), without the Statistics Toolbox.
    x = sort(x(:));  n = numel(x);
    if n == 0, y = NaN; return; end
    pos = 100 * ((1:n)' - 0.5) / n;
    if p <= pos(1), y = x(1); elseif p >= pos(end), y = x(end);
    else, y = interp1(pos, x, p); end
end
% Author: Diellor Basha, 2026
