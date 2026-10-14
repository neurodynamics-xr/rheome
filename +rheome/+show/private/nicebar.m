function nicebar(cb, ax)
% NICEBAR  A one-sided colourbar on round limits with 3-6 round ticks (tokens: axis.ticks), both ends
% labelled. The range is only widened, so no data leaves the scale. Diverging bars keep -L, 0, +L.
%
% Author: Diellor Basha, 2026
    lim = ax.CLim;  r = diff(lim);
    if ~(r > 0), return; end
    d = 10^floor(log10(r));
    for s = [0.2 0.5 1 2 5] * d                  % the finest round step giving at most 4 ticks (a key-band bar is short)
        if floor(r/s) + 1 <= 4, break; end
    end
    k = floor(lim(1)/s + 1e-9):ceil(lim(2)/s - 1e-9);
    tk = k * s;  clim(ax, tk([1 end]));
    cb.Ticks = tk;  cb.TickLabels = compose('%g', tk);  cb.TickLabelsMode = 'manual';
    p = strsplit(cb.Tag, '|');                    % keep the house tag (read by figure_check) in step
    if numel(p) >= 5, p{3} = sprintf('%g', tk(1));  p{4} = sprintf('%g', tk(end));  cb.Tag = strjoin(p, '|'); end
end
