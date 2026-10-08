function drawScalogram(app)
% DRAWSCALOGRAM  Spatial scale against time for the whole page.
%
% ⭐ THIS IS WHERE THE HYPOTHESIS IS VISIBLE. If coherent patterns fragment at
% desynchronisation, energy MIGRATES UP the scale axis -- a shape you can see, not a slope
% you have to trust. Band power sits directly beneath on the same time axis, so scale and
% power are read together.
%
% ⚠ COLUMN-NORMALISED, DELIBERATELY. Each frame's scales are shown as FRACTIONS summing to
% one, so the image shows how energy is DISTRIBUTED across scales independently of how much
% there is. Without it a burst dominates every column and the picture becomes a restatement
% of the power trace -- the exact confusion the analysis exists to avoid. Absolute power is
% not lost; it is the panel below.
%
% ⚠ THE BANK'S MEMBER ORDER IS NOT MONOTONIC IN SCALE. Member 1 is the LOWPASS (127 mm here)
% and the wavelets then run FINE to COARSE (36, 39, 45, 59, 79, 99 mm). Plotting by member
% index would scramble the y axis into something that looks like structure and is pure
% ordering -- the same trap as the unsorted lambda axis. Rows are sorted by wavenumber:
% coarse at the bottom, fine at the top, so "fragmenting" reads as "moving up".
%
% ⚠ imagesc ASSUMES A UNIFORM AXIS, and the scales are log-spaced. Rows are drawn by INDEX
% and the wavelengths go in the tick LABELS.
%
% Author: Diellor Basha, 2026

    fp = app.Page;
    ax = app.ScalogramAx;  cla(ax);

    E  = app.ScalogramData;                                % [nScale x nT], raw
    kc = double(reshape(centerWavenumbers(fp.GraphBank), [], 1));
    [kS, ord] = sort(kc, 'ascend');                        % small k = coarse = bottom row
    E  = E(ord, :);
    Efrac = E ./ max(sum(E, 1), realmin);

    imagesc(ax, [fp.Time(1) fp.Time(end)], [1 fp.NumScales], Efrac);
    set(ax, 'YDir', 'normal');
    colormap(ax, parula(256));
    ax.YTick = 1:fp.NumScales;
    ax.YTickLabel = compose('%.0f', 1000*2*pi./kS);
    ylabel(ax, 'scale (mm)   coarse \rightarrow fine \uparrow');
    xlim(ax, [fp.Time(1) fp.Time(end)]);

    % the centroid in the SAME row coordinates -- the summary drawn on the distribution it
    % summarises, so it is obvious when one number is failing to represent the other
    rowC = (1:fp.NumScales) * Efrac;
    hold(ax, 'on');
    plot(ax, fp.Time, rowC, '-', 'Color', [1 1 1 0.8], 'LineWidth', 1.2);
    hold(ax, 'off');
    app.shadeMargin(ax);
    title(ax, 'cortical scalogram — fraction of energy per scale   (white: centroid row)');
end

% Author: Diellor Basha, 2026
