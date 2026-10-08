function out = sensor_limits(outDir)
% DEMOS.SENSOR_LIMITS  Which wavelengths does each array admit? Predicted, then verified.
%
%   rheome.demos.sensor_limits()
%   out = rheome.demos.sensor_limits('_figures')
%
%   1 windows   the admitted band per array, in millimetres of sensor separation
%   2 sweep     planted wavelength against recovered, until measurement fails
%   3 oblique   a wave crossing a probe at angle theta, read as k*cos(theta)
%   4 topology  which arrays carry a winding number at all
%
% ⭐ THE WINDOW IS PREDICTED BEFORE ANY SIGNAL EXISTS. Figure 1 comes from coordinates
% alone: 2*pitch at the bottom, aperture at the top. Figure 2 then plants waves and sweeps
% across both edges until recovery breaks, so the bars in figure 1 are TESTED rather than
% drawn.
%
% ⭐ WAVELENGTH IS RECOVERED SPECTRALLY, WHICH IS WHY THE PROBE APPEARS AT ALL. The route
% here is lambda_hat -> k = sqrt(lambda_hat/alpha), which needs the operator and nothing
% else -- no faces, no triangulation. So the probe is measured by the same code as the
% grids, and where it differs it differs for a stated reason rather than a missing branch.
%
% ⭐ THE PROBE'S TWO LIMITS ARE DIFFERENT IN KIND, and figures 3 and 4 separate them.
% Figure 3: a wave crossing at theta projects, so the probe reads k*cos(theta) -- a longer
% apparent wavelength and a faster apparent speed, biased and uncorrectable from the probe
% alone. Figure 4: a winding number is summed around a face (a 2-cell), and a chain admits
% no faces, so it is not a hard measurement on it, it is an undefined one.
%
% ⚠ EVERY LENGTH IS SENSOR SEPARATION, and nothing here is a claim about sources.
%
% INPUTS:
%   outDir  '' for screen, or a folder for PNGs
% OUTPUT:
%   out  .ok .checks .arrays .win .sweep .obliqueK .obliqueTheta
%
% See also: rheome.demos.sensor_graph, rheome.sensors.window, rheome.sensors.calibrate
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});

    % ⚠ 64 SITES, NOT 32. Figure 3 plants a wave and projects it by cos(theta); at theta =
    % 60 degrees the apparent wavelength DOUBLES. On a 32-site, 100 um-pitch probe the
    % aperture is 31*100um = 3.1 mm, and doubling the figure-3 planted wavelength (see
    % below) would push it past that aperture -- the test would then be measuring the
    % array's ceiling, not the cos(theta) projection it claims to measure. 64 sites gives a
    % 63*100um = 6.3 mm aperture, comfortably clearing it.
    arrays = {rheome.sensors.linear('NumSites', 64, 'Pitch', 100e-6), ...
              rheome.sensors.grid('utah'), ...
              rheome.sensors.grid('ecog'), ...
              rheome.sensors.eeg('NumChannels', 64)};
    nA  = numel(arrays);
    win = cellfun(@rheome.sensors.window, arrays, 'UniformOutput', false);

    % ⭐ ONE OPERATOR AND ONE CALIBRATION PER ARRAY, SHARED ACROSS FIGURES. Figures 1-3 all
    % need rheome.sensors.graph/sensors.calibrate on the same arrays; hoisting them here (as
    % rheome.demos.sensor_graph already does) means figure 1's bar and figure 2's sweep read the
    % SAME calibrated object rather than two independently-recomputed ones that would only
    % coincidentally describe the same array.
    graphs = cellfun(@rheome.sensors.graph, arrays, 'UniformOutput', false);
    cal    = cellfun(@rheome.sensors.calibrate, graphs, 'UniformOutput', false);

    fprintf('\n=== rheome.demos.sensor_limits : %d arrays ===\n', nA);

    % ---------- Fig 1: the predicted windows ----------
    f1 = fbd_fig(doExport, [60 60 900 420]);
    hold on;
    for i = 1:nA
        cw = cal{i};
        lo = win{i}.LambdaMin*1e3;  us = cw.lambdaUsable*1e3;  hi = win{i}.LambdaMax*1e3;
        plot([lo hi], [i i], '-', 'LineWidth', 3, 'Color', [.75 .75 .75]);
        if isnan(us)
            % No usable band. Say so on the axis rather than drawing a bar that is not there.
            text(hi*1.15, i, sprintf('NO usable band (%.1f pitches across)', ...
                arrays{i}.Aperture/arrays{i}.Pitch), 'FontSize', 9, 'Color', [.6 0 0]);
        else
            plot([us hi], [i i], '-', 'LineWidth', 9);
            text(hi*1.15, i, sprintf('usable %.2g - %.3g mm', us, hi), 'FontSize', 9);
        end
    end
    set(gca, 'XScale', 'log', 'YTick', 1:nA, ...
        'YTickLabel', cellfun(@(a) a.Name, arrays, 'UniformOutput', false), ...
        'TickLabelInterpreter', 'none');
    xlim([0.1 1e3]); ylim([0.4 nA+0.6]); grid on;
    xlabel('wavelength (mm of sensor separation)');
    title(['the admitted band. Grey = down to the alias floor; thick = where the ' ...
           'quadratic model still holds']);
    fbd_save(f1, outDir, 'sensor_limits_1_windows');

    % ---------- Fig 2: sweep until it breaks ----------
    sweep = struct('name', {}, 'lamPlanted', {}, 'lamRecovered', {}, 'inWindow', {});
    f2 = fbd_fig(doExport, [70 70 1200 420]);
    for i = 1:nA
        arr = arrays{i};
        G   = graphs{i};
        c   = cal{i};
        dir1 = c.Directions(1, :);

        lamP = logspace(log10(0.5*win{i}.LambdaMin), log10(3*win{i}.LambdaMax), 30);
        lamR = nan(size(lamP));
        for j = 1:numel(lamP)
            kv = (2*pi/lamP(j)) * dir1;
            u  = rheome.sensors.sample(arr, @(P) cos(P * kv.'));
            u  = u - mean(u);
            nu = u.' * u;
            % ⚠ AT THE ALIAS FLOOR u GOES CONSTANT. A wave of exactly one pitch per cycle
            % samples the same phase at every sensor, so u - mean(u) is identically zero and
            % there is no wavenumber left to read. NaN is the honest answer, not a number.
            if nu <= 1e-12 * arr.nCh, continue; end
            lh = (u.' * (G.L * u)) / nu;
            lamR(j) = 2*pi / sqrt(max(lh, eps) / c.alpha);
        end
        % ⭐ THE USABLE WINDOW, NOT THE ALIAS WINDOW. lambda ~ alpha k^2 is a small-k
        % expansion and is already 10% low several sensors per cycle above the alias floor,
        % so c.lambdaUsable -- not 2*pitch -- is where accurate recovery actually starts.
        % ⚠ lambdaUsable IS NaN WHEN THE ARRAY HAS NO USABLE BAND AT ALL, and NaN in a
        % comparison yields all-false -- an empty band that a later max() turns into -Inf,
        % passing vacuously. Record the emptiness explicitly: an array with no usable band
        % is a RESULT this demo reports, not a case to skip quietly.
        if isnan(c.lambdaUsable)
            inWin = false(size(lamP));
        else
            inWin = lamP >= c.lambdaUsable & lamP <= win{i}.LambdaMax;
        end
        sweep(i) = struct('name', arr.Name, 'lamPlanted', lamP, ...
                          'lamRecovered', lamR, 'inWindow', inWin); %#ok<AGROW>

        subplot(1, nA, i);
        loglog(lamP*1e3, lamR*1e3, '.-'); hold on;
        loglog(lamP*1e3, lamP*1e3, 'k--');
        xline(win{i}.LambdaMin*1e3, 'r--');            % alias floor
        if ~isnan(c.lambdaUsable)
            xline(c.lambdaUsable*1e3, 'r-');           % where the model is still good
        end
        xline(win{i}.LambdaMax*1e3, 'r-');             % aperture ceiling
        grid on; axis tight;
        xlabel('planted \lambda (mm)'); ylabel('recovered \lambda (mm)');
        title(arr.Name, 'Interpreter', 'none');
    end
    sgtitle(['sensor\_limits -- Fig 2: the bars in Fig 1, tested. Solid red = the USABLE ' ...
             'window, dashed = the alias floor below it']);
    fbd_save(f2, outDir, 'sensor_limits_2_sweep');

    % ---------- Fig 3: the oblique wave on a probe ----------
    probe = arrays{1};
    Gp    = graphs{1};
    cp    = cal{1};
    axisDir = cp.Directions(1, :);
    perp    = null(axisDir).';  perp = perp(1, :);
    theta   = linspace(0, pi/3, 7);
    % ⚠ 1.5 mm. R3 (64 sites, not 32) is what makes this figure work at all: on the 32-site
    % probe R3 replaced, EITHER planted wavelength fails the 5% self-check below (1mm:
    % 5.75%, 1.5mm: 6.94%, both measured, not estimated). On the 64-site probe both 1mm and
    % 1.5mm pass; the choice between them is a marginal improvement, not a fix. On the
    % test's own RelTol metric (which divides by cos(theta), so it is STRICTER at theta =
    % 60 where cos is smallest), 1.5mm measures 2.10% against 1mm's 2.22%, and 1.5mm's
    % sweep to 3.0mm (at theta = 60, lambda/cos(60) doubles it) has geometric mean 2.12mm --
    % close to the geometric centre of this probe's own usable band, sqrt(0.75mm * 6.3mm)
    % = 2.18mm. lam0 = 2mm, by contrast, measures 5.58% and fails outright.
    lam0    = 1.5e-3;
    kmag    = 2*pi/lam0;
    kObs    = zeros(size(theta));
    for j = 1:numel(theta)
        kv = kmag * (cos(theta(j))*axisDir + sin(theta(j))*perp);
        u  = rheome.sensors.sample(probe, @(P) cos(P * kv.'));
        u  = u - mean(u);
        lh = (u.' * (Gp.L * u)) / (u.' * u);
        kObs(j) = sqrt(max(lh, eps) / cp.alpha);
    end

    f3 = fbd_fig(doExport, [80 80 1000 400]);
    subplot(1, 2, 1);
    plot(theta*180/pi, kObs/kObs(1), 'o-'); hold on;
    plot(theta*180/pi, cos(theta), 'k--'); grid on;
    xlabel('\theta from the shank axis (deg)'); ylabel('k_{observed} / k_{true}');
    legend({'measured', 'cos\theta'}, 'Location', 'southwest');
    title('a probe sees only the projection');
    subplot(1, 2, 2);
    plot(theta*180/pi, 1./cos(theta), 'o-'); grid on;
    xlabel('\theta (deg)'); ylabel('apparent speed / true');
    title('so any speed it reports is biased HIGH');
    sgtitle('sensor\_limits -- Fig 3: c/cos\theta, and nothing on the probe can correct it');
    fbd_save(f3, outDir, 'sensor_limits_3_oblique');

    % ---------- Fig 4: which arrays carry a winding number ----------
    f4 = fbd_fig(doExport, [90 90 900 380]);
    hasFac = cellfun(@(w) w.HasFaces, win);
    bar(double(hasFac)); grid on; ylim([0 1.3]);
    set(gca, 'XTickLabel', cellfun(@(a) a.Name, arrays, 'UniformOutput', false), ...
        'TickLabelInterpreter', 'none', 'YTick', [0 1], ...
        'YTickLabel', {'undefined', 'defined'});
    xtickangle(20);
    title('winding number: summed around a face, and a chain admits none to wind around');
    fbd_save(f4, outDir, 'sensor_limits_4_topology');

    % ---------- self-checks ----------
    fprintf('\n  self-checks\n');
    for i = 1:nA
        s = sweep(i);
        % ⚠ GUARD THE EMPTY BAND. When c.lambdaUsable is NaN, s.inWindow is all-false and
        % max() on an empty selection returns [] -- fbd_check then evaluates
        % isfinite([]) && ..., which is not a scalar logical and THROWS
        % (MATLAB:nonLogicalConditional). Loud rather than silent, but this function twice
        % takes trouble to handle the no-usable-band case gracefully elsewhere; do the same
        % here. Recording NaN, not skipping, keeps the missing band a FAILING check (per
        % fbd_check's own "NaN is a failure" contract) rather than one that vanishes.
        if any(s.inWindow)
            err = max(abs(s.lamRecovered(s.inWindow)./s.lamPlanted(s.inWindow) - 1));
        else
            err = NaN;
        end
        chk(end+1) = fbd_check(sprintf('%s: in-window recovery', s.name), ...
            err, 0, 0.15, 'rel'); %#ok<AGROW>
    end
    chk(end+1) = fbd_check('probe: cos(theta) projection', ...
        max(abs(kObs./kObs(1) - cos(theta))), 0, 0.05, '');
    % ⚠ NOT win{1}.HasFaces. window.m defines hasFaces = arr.Dim >= 2, so checking it here
    % would only restate that a Dim-1 array has Dim < 2 -- a tautology that cannot catch a
    % regression in rheome.sensors.graph, the component that could actually triangulate faces onto
    % a chain. Checking Gp.Faces (Gp = graphs{1}, the probe's own constructed operator) instead
    % measures the graph itself: rheome.sensors.graph triangulates faces only for Dim >= 2 arrays
    % (tSensorGeometry: "aProbeHasNoFacesAndThatIsTheResultNotAGap"), so a future regression
    % that accidentally triangulated a 1D array would show up here as a nonzero face count,
    % where the old check could not have noticed.
    chk(end+1) = fbd_check('probe: no readout faces (winding undefined)', ...
        size(Gp.Faces, 1), 0, 0.5, '');
    [ok, ~] = fbd_report(chk);

    out = struct('ok', ok, 'checks', chk, 'arrays', {arrays}, 'win', {win}, ...
                 'sweep', sweep, 'obliqueK', kObs, 'obliqueTheta', theta);
end

% Author: Diellor Basha, 2026
