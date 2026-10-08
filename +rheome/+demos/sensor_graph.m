function out = sensor_graph(outDir)
% DEMOS.SENSOR_GRAPH  From sensor coordinates to a calibrated operator, in five figures.
%
%   rheome.demos.sensor_graph()
%   out = rheome.demos.sensor_graph('_figures')     % also export the figures as PNGs
%
%   1 geometry     five arrays: probe, Utah, ECoG, MEG helmet, EEG cap
%   2 graph        weights, degrees, and the readout triangulation
%   3 spectrum     lambda and the low eigenmodes
%   4 calibration  alpha by two routes, and the departure that is the resolution floor
%   5 routes       Chebyshev against the exact eigenbasis on the same operator
%
% ⭐ WHAT THIS DEMO EARNS. A Gaussian-kernel graph Laplacian has dimensionless lambda, so on
% its own it cannot report a wavelength in millimetres or a speed in metres per second.
% Figure 4 recovers alpha in lambda ~ alpha k^2 by two independent routes -- a closed form
% and a regression -- and their agreement is what puts units on everything downstream.
%
% ⭐ THE DEPARTURE IS THE MEASUREMENT. In figure 4, lambda_hat(k) leaves alpha*k^2 as the
% wavelength approaches twice the pitch. That is not error: it is the array's resolution
% floor, read off a curve rather than quoted from a rule of thumb.
%
% ⚠ EVERYTHING HERE IS SENSOR SPACE. Lengths are millimetres of separation between sensors.
% Nothing in this demo is a claim about what produced the signal.
%
% INPUTS:
%   outDir  '' to show figures on screen, or a folder to export PNGs
% OUTPUT:
%   out  .ok .checks .arrays .cal .alphaAgree .regular .hasWindow .handled .chebErr
%
% See also: rheome.demos.sensor_limits, rheome.sensors.calibrate, rheome.sensors.window
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});

    arrays = {rheome.sensors.linear('NumSites', 64, 'Pitch', 100e-6), ...
              rheome.sensors.grid('utah'), ...
              rheome.sensors.grid('ecog'), ...
              i_meg(), ...
              rheome.sensors.eeg('NumChannels', 64)};
    arrays = arrays(~cellfun(@isempty, arrays));
    nA     = numel(arrays);
    regular = cellfun(@(a) any(strcmp(a.Kind, {'linear','grid'})), arrays);

    fprintf('\n=== rheome.demos.sensor_graph : %d arrays ===\n', nA);

    graphs = cell(1, nA);  cal = cell(1, nA);
    for i = 1:nA
        graphs{i} = rheome.sensors.graph(arrays{i});
        cal{i}    = rheome.sensors.calibrate(graphs{i});
    end
    alphaAgree = cellfun(@(c) c.agree, cal);
    hasWindow  = cellfun(@(c) c.hasFitWindow, cal);

    % ---------- Fig 1: the geometries ----------
    f1 = fbd_fig(doExport, [60 60 1400 340]);
    for i = 1:nA
        subplot(1, nA, i);
        P = arrays{i}.Pos * 1e3;
        plot3(P(:,1), P(:,2), P(:,3), '.', 'MarkerSize', 9); axis equal; grid on;
        if arrays{i}.Dim == 1, view(0, 0); else, view(3); end
        title(sprintf('%s\n%d ch, pitch %.3g mm', arrays{i}.Name, arrays{i}.nCh, ...
                      arrays{i}.Pitch*1e3), 'Interpreter', 'none');
        xlabel('mm');
    end
    sgtitle('sensor\_graph -- Fig 1: five arrays, coordinates only');
    fbd_save(f1, outDir, 'sensor_graph_1_geometry');

    % ---------- Fig 2: the graphs ----------
    f2 = fbd_fig(doExport, [70 70 1400 400]);
    for i = 1:nA
        subplot(2, nA, i);
        spy(graphs{i}.W, 3);
        title(sprintf('W  (K = %d)', graphs{i}.K));
        subplot(2, nA, nA + i);
        histogram(full(sum(graphs{i}.W, 2)), 16); grid on;
        xlabel('degree');
        if isempty(graphs{i}.Faces)
            title('no faces -- winding undefined');
        else
            title(sprintf('%d readout faces', size(graphs{i}.Faces, 1)));
        end
    end
    sgtitle('sensor\_graph -- Fig 2: Gaussian-weighted kNN, and the readout loops');
    fbd_save(f2, outDir, 'sensor_graph_2_graph');

    % ---------- Fig 3: the spectra ----------
    modesAll = cell(1, nA);
    f3 = fbd_fig(doExport, [80 80 1400 400]);
    subplot(1, 4, 1); hold on;
    for i = 1:nA
        modesAll{i} = rheome.sensors.modes(graphs{i});
        plot(modesAll{i}.Lambda / modesAll{i}.Lmax, '.-');
    end
    grid on; xlabel('mode index'); ylabel('\lambda / \lambda_{max}');
    legend(cellfun(@(a) a.Name, arrays, 'UniformOutput', false), 'Location', 'southeast', ...
           'Interpreter', 'none');
    title('normalised spectra');
    iG = find(cellfun(@(a) strcmp(a.Kind,'grid'), arrays), 1, 'last');
    for j = 1:3
        subplot(1, 4, j + 1);
        P = arrays{iG}.Pos * 1e3;
        scatter(P(:,1), P(:,2), 26, modesAll{iG}.Phi(:, j+1), 'filled');
        axis equal off; title(sprintf('%s: mode %d', arrays{iG}.Name, j+1), ...
                              'Interpreter', 'none');
    end
    sgtitle('sensor\_graph -- Fig 3: the spectrum, and what its low modes look like');
    fbd_save(f3, outDir, 'sensor_graph_3_spectrum');

    % ---------- Fig 4: the calibration ----------
    f4 = fbd_fig(doExport, [90 90 1400 420]);
    subplot(1, 3, 1); hold on;
    for i = 1:nA
        loglog(cal{i}.kGrid, cal{i}.lambdaHat, '.-');
    end
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
    xlabel('k (rad/m)'); ylabel('$\hat\lambda$', 'Interpreter', 'latex');
    title('Rayleigh quotient vs planted k');
    legend(cellfun(@(a) a.Name, arrays, 'UniformOutput', false), 'Location', 'northwest', ...
           'Interpreter', 'none');

    % ⚠ LOG y-AXIS. alpha spans five decades across these arrays (1e-8 to 1.6e-3 m^2) --
    % on a linear axis the laminar probe's bars are 0.0006% of full scale, i.e. invisible.
    % The two routes' AGREEMENT (in %) is annotated above each bar pair because that is
    % the number this panel exists to show, and a bar-height comparison alone cannot
    % convey it once the axis spans five decades.
    subplot(1, 3, 2); hold on;
    alphaReg = cellfun(@(c) c.alpha, cal);
    alphaCl  = cellfun(@(c) c.alphaClosed, cal);
    bar([alphaReg; alphaCl].');
    set(gca, 'YScale', 'log'); grid on;
    % A newline inside a rotated XTickLabel renders unreliably under software OpenGL
    % (the second line drifts onto the neighbouring tick) -- one line per label instead.
    tickLabels = cell(1, nA);
    for i = 1:nA
        tickLabels{i} = sprintf('%s (%s, nFit=%d)', arrays{i}.Name, ...
            fbd_ternary(hasWindow(i), 'window', 'NO window'), cal{i}.nFit);
    end
    set(gca, 'XTick', 1:nA, 'XTickLabel', tickLabels, 'TickLabelInterpreter', 'none', ...
        'FontSize', 8);
    xtickangle(30);
    legend({'regression', 'closed form'}, 'Location', 'best');
    for i = 1:nA
        text(i, 1.25 * max(alphaReg(i), alphaCl(i)), sprintf('%.0f%%', 100*alphaAgree(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    ylabel('\alpha (m^2)'); title('two routes to \alpha  (label = %agree)');

    % ⚠ THE xline MARKS lambdaUsable, NOT AN EYEBALLED KNEE. It is where the quadratic
    % model itself stops tracking lambda_hat to 10% (rheome.sensors.calibrate), so the reader
    % sees where the resolution floor actually sits rather than guessing it off the curve.
    % NaN means no usable band exists on that array at all -- skip the line and say so.
    subplot(1, 3, 3); hold on;
    ph = gobjects(1, nA);
    for i = 1:nA
        lamN = 2*pi ./ cal{i}.kGrid / arrays{i}.Pitch;          % sensors per cycle
        ph(i) = plot(lamN, cal{i}.lambdaHat ./ (cal{i}.alpha * cal{i}.kGrid.^2), '.-');
    end
    set(gca, 'XScale', 'log'); grid on; yline(1, 'k--', 'HandleVisibility', 'off');
    for i = 1:nA
        if isfinite(cal{i}.kUsable)
            xline(cal{i}.lambdaUsable / arrays{i}.Pitch, '--', 'Color', ph(i).Color, ...
                'HandleVisibility', 'off');
        else
            fprintf('  %-38s lambdaUsable is NaN -- no usable band on this array\n', ...
                arrays{i}.Name);
        end
    end
    legend(ph, cellfun(@(a) a.Name, arrays, 'UniformOutput', false), 'Location', 'best', ...
           'Interpreter', 'none');
    xlabel('sensors per wavelength'); ylabel('$\hat\lambda\,/\,\alpha k^2$', 'Interpreter', 'latex');
    title('the departure IS the resolution floor  (dashed = \lambda_{usable})');
    sgtitle('sensor\_graph -- Fig 4: \lambda carries metres once \alpha is measured');
    fbd_save(f4, outDir, 'sensor_graph_4_calibration');

    % ---------- Fig 5: the two routes on one operator ----------
    Gc  = graphs{iG};
    Me  = rheome.sensors.modes(Gc, 'Route', 'eigen');
    Mc  = rheome.sensors.modes(Gc, 'Route', 'chebyshev', 'Order', 80);
    g   = @(l) exp(-0.4 * l / Me.Lmax);
    X   = randn(Gc.nV, 1);
    Ye  = Me.Phi * (g(Me.Lambda) .* (Me.Phi.' * X));
    Yc  = Mc.T.filter(g, X);
    chebErr = norm(Yc - Ye) / norm(Ye);

    f5 = fbd_fig(doExport, [100 100 1100 400]);
    subplot(1, 2, 1); plot(Ye, Yc, '.'); grid on; axis equal; hold on;
    lims = [min([Ye; Yc]) max([Ye; Yc])];
    plot(lims, lims, 'k--');                  % NOT refline -- that is Statistics Toolbox
    xlabel('eigen route'); ylabel('Chebyshev route');
    title('the same operator, two ways');
    subplot(1, 2, 2);
    P = arrays{iG}.Pos * 1e3;
    scatter(P(:,1), P(:,2), 30, Yc - Ye, 'filled'); axis equal off; colorbar;
    title(sprintf('difference, rel. err %.2g', chebErr));
    sgtitle('sensor\_graph -- Fig 5: Chebyshev is the route that SCALES, and it agrees');
    fbd_save(f5, outDir, 'sensor_graph_5_routes');

    % ---------- self-checks ----------
    fprintf('\n  self-checks\n');
    % handled is populated BY the loops below, at the point each array is actually dealt
    % with -- NOT recomputed from regular/hasWindow afterwards. A test that recomputed the
    % loop conditions would be a tautology over two booleans and would pass even if a loop
    % were deleted; this vector only becomes true where a loop body actually ran.
    handled = false(1, nA);
    % ⚠ 0.20, NOT 0.05. The closed form assumes every sensor sees the same neighbourhood,
    % and these arrays are SMALL -- a Utah array is 36% boundary and an ECoG grid 44%. The
    % looser bound here against the 0.05 that tSensorCalibrate holds a 32x32 lattice to is
    % the boundary effect, and Fig 4 shows it rather than hiding it.
    for i = find(regular & hasWindow)
        chk(end+1) = fbd_check(sprintf('%s: alpha routes agree', arrays{i}.Name), ...
            cal{i}.agree, 0, 0.20, 'rel'); %#ok<AGROW>
        handled(i) = true;
    end
    % Arrays WITHOUT the window are reported, not scored -- the absence IS the finding.
    for i = find(~hasWindow)
        fprintf('  %-38s NO small-k fit window (%.1f pitches across, nFit = %d)\n', ...
            arrays{i}.Name, arrays{i}.Aperture/arrays{i}.Pitch, cal{i}.nFit);
        handled(i) = true;
    end
    % Irregular arrays WITH a window are reported too -- the disagreement IS the error bar.
    % The closed form assumes every sensor sees the same neighbourhood; a helmet or cap does
    % not, so a real fit window does not make the two routes trustworthy against each other
    % the way it does on a lattice. Not scored against the 0.20 lattice bound -- that bound
    % is about the CODE being right on a geometry where the closed form is exact in the
    % interior, and an irregular array's disagreement is the geometry's fact, not the code's.
    for i = find(~regular & hasWindow)
        fprintf(['  %-38s alpha routes DISAGREE by %.0f%% (fit window OK, nFit = %d) -- ' ...
                 '%.0f%% on any speed from this array, alphaSpread = %.2f\n'], ...
            arrays{i}.Name, 100*cal{i}.agree, cal{i}.nFit, ...
            100*(sqrt(cal{i}.alphaClosed/cal{i}.alpha) - 1), cal{i}.alphaSpread);
        handled(i) = true;
    end
    chk(end+1) = fbd_check('Chebyshev vs eigen', chebErr, 0, 1e-3, 'rel');
    for i = 1:nA
        chk(end+1) = fbd_check(sprintf('%s: lambda_1 = 0', arrays{i}.Name), ...
            modesAll{i}.Lambda(1), 0, 1e-8, ''); %#ok<AGROW>
    end
    [ok, ~] = fbd_report(chk);

    out = struct('ok', ok, 'checks', chk, 'arrays', {arrays}, 'cal', {cal}, ...
                 'alphaAgree', alphaAgree, 'regular', regular, 'hasWindow', hasWindow, ...
                 'handled', handled, 'chebErr', chebErr);
end

function arr = i_meg()
% The one array with real coordinates -- skipped cleanly when the cache is absent, exactly
% as rheome.demos.filterbank_cwt skips without the Wavelet Toolbox.
    arr = [];
    try
        arr = rheome.sensors.meg('subject01');
    catch
        fprintf('  (subject01 not cached -- the MEG array is skipped)\n');
    end
end

% Author: Diellor Basha, 2026
