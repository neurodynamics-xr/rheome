function plainlog(ax)
% PLAINLOG  Plain-number tick labels on the log axes of ax: MATLAB writes 10^k with the exponent at 0.8x
% the tick size, below the print floor of the figure tokens (6 pt).
%
% Author: Diellor Basha, 2026
    for a = 'XY'
        if strcmp(ax.([a 'Scale']), 'log'), ax.([a 'TickLabel']) = compose('%g', ax.([a 'Tick'])); end
    end
end
