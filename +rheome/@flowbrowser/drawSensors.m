function drawSensors(app)
% DRAWSENSORS  A few band-filtered sensor traces under the cortex, with a time cursor.
%
% ⭐ WHY THE BAND-FILTERED SIGNAL, NOT THE RAW. These are Re(SensorCoefficients) -- the
% band's contribution at each sensor, i.e. exactly what the flow kernel is being applied to.
% Raw traces would be dominated by out-of-band activity and would not correspond to what the
% surfaces show.
%
% Channels are the highest-variance ones over the page's core, so the traces carry the band
% rather than showing whichever sensors happen to come first in the array.
%
% ⚠ TRACES ARE DRAWN ONCE PER PAGE. Only the cursor moves per frame -- see updateCursors.
%
% Author: Diellor Basha, 2026

    fp = app.Page;
    ax = app.SensorAx;  cla(ax);

    X = real(double(fp.SensorCoefficients));            % [C x nT]
    v = var(X(:, fp.Core), 0, 2);
    [~, ord] = sort(v, 'descend');
    nShow = min(app.NumSensorTraces, size(X,1));
    pick  = sort(ord(1:nShow));

    % stack them on a common offset so shape is comparable across channels
    s = median(std(X(pick, fp.Core), 0, 2));
    if ~(s > 0), s = 1; end
    hold(ax, 'on');
    for k = 1:nShow
        plot(ax, fp.Time, X(pick(k), :)/(4*s) + k, 'LineWidth', 0.6);
    end
    hold(ax, 'off');

    ax.YTick = 1:nShow;
    ax.YTickLabel = app.Bundle.pager.ChannelName(app.Bundle.chSel(pick));
    ylim(ax, [0.2 nShow + 0.8]);
    xlim(ax, [fp.Time(1) fp.Time(end)]);
    app.shadeMargin(ax);
    grid(ax, 'on');
    xlabel(ax, 'time (s)');
    title(ax, sprintf('%d strongest sensors, %s band (%.1f–%.1f Hz)', ...
        nShow, app.BandName, min(fp.CenterFrequencies), max(fp.CenterFrequencies)));

    app.SensorCursor = xline(ax, fp.Time(app.Frame), 'r-', 'LineWidth', 1.2);
end

% Author: Diellor Basha, 2026
