function out = sensor_dynamics(outDir)
% DEMOS.SENSOR_DYNAMICS  What a moving pattern looks like on an array, and in its spectrum.
%
%   rheome.demos.sensor_dynamics()
%   out = rheome.demos.sensor_dynamics('_figures')     % also export PNGs and the web JSON
%
%   1 structure  one travelling wave: snapshots, sensor traces, carpet, joint spectrum
%   2 instruments the SAME wave through four geometries -- what each one sees
%   3 families   travelling / rotating / standing / diffusing, side by side
%   4 the claim  travelling vs standing: identical per-channel, opposite in the spectrum
%   5 dispersion a wave train of one speed -> a ridge -> the speed back out
%
% ⭐ WHY FIGURE 4 IS THE ARGUMENT. A travelling wave and a standing wave put their energy at
% the SAME (lambda, omega) point -- same spatial scale, same frequency -- and produce nearly
% the same trace at any single sensor. Nothing in a per-channel analysis separates them. What
% separates them is the COMPLEX SPATIAL PATTERN at that frequency: standing collapses every
% sensor onto one line through the origin (a common phase), travelling spreads them around a
% circle (a phase that advances with position). That is a property of the field across the
% array, which is exactly what a graph-spectral method sees and a channel-by-channel one
% cannot.
%
% ⚠ EVERY WAVE IS PLANTED AT A WAVELENGTH THAT ARRAY CAN ACTUALLY SEE. Each geometry gets a
% wave at the geometric centre of its own usable band (rheome.sensors.calibrate). Planting one
% physical wavelength across all four would be measuring the instruments' windows, which
% rheome.demos.sensor_limits already does, rather than their dynamics.
%
% ⚠ SENSOR SPACE ONLY. Speeds are metres per second ACROSS THE ARRAY and wavelengths are
% millimetres of SENSOR SEPARATION. Nothing here is a claim about what produced a signal.
%
% INPUTS:
%   outDir  '' for screen, or a folder for PNGs + sensor_dynamics.json
% OUTPUT:
%   out  .ok .checks .structure .families .standingIndex .speedFit
%
% See also: rheome.demos.sensor_graph, rheome.demos.sensor_limits, rheome.sensors.calibrate
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});
    fs = 200;  nT = 256;  tv = (0:nT-1)/fs;      % 1.28 s at 200 Hz

    fprintf('\n=== rheome.demos.sensor_dynamics : %d samples at %g Hz ===\n', nT, fs);

    % ---------- the instruments ----------
    arrays = {rheome.sensors.linear('NumSites',64,'Pitch',100e-6), ...
              rheome.sensors.grid('utah'), rheome.sensors.grid('ecog'), i_meg()};
    arrays = arrays(~cellfun(@isempty, arrays));
    nA = numel(arrays);
    G = cell(1,nA); cal = cell(1,nA); pl = cell(1,nA); lamP = zeros(1,nA);
    for i = 1:nA
        G{i}   = rheome.sensors.graph(arrays{i});
        cal{i} = rheome.sensors.calibrate(G{i});
        pl{i}  = sen_plane(arrays{i});
        lamP(i) = i_pickwavelength(arrays{i}, cal{i});
    end
    iG = find(cellfun(@(a) strcmp(a.Kind,'grid'), arrays), 1, 'last');   % the ECoG grid

    % ---------- Fig 1: the make-up of one travelling wave ----------
    f0 = 10;                                        % Hz
    [U1, tru1] = i_travelling(arrays{iG}, pl{iG}, lamP(iG), f0, tv);
    an = i_analyse(U1, G{iG}, cal{iG}, fs);

    f1 = fbd_fig(doExport, [50 50 1400 780]);
    snapIdx = round(linspace(1, round(fs/f0), 5));  % one temporal period
    for j = 1:5
        subplot(3,5,j);
        i_scatterfield(arrays{iG}, U1(:,snapIdx(j)), max(abs(U1(:))));
        title(sprintf('t = %.0f ms', 1e3*tv(snapIdx(j))), 'FontWeight','normal');
    end
    subplot(3,1,2);
    i_traces(U1, pl{iG}, arrays{iG}, tv, 6);
    title('sensor traces along the propagation direction -- the same wave, later at each');
    subplot(3,2,5);
    i_carpet(U1, pl{iG}, arrays{iG}, tru1.khat, tv);
    title('carpet: sensors ordered along k, time across');
    subplot(3,2,6);
    i_jointimage(an);
    title(sprintf('joint spectrum: one point at (k, f) = (%.0f rad/m, %.1f Hz)', ...
        an.kPeak, an.fPeak));
    sgtitle(sprintf(['sensor\\_dynamics -- Fig 1: a %.0f mm wave at %.0f Hz on the %s ' ...
        '(planted speed %.2f m/s)'], lamP(iG)*1e3, f0, arrays{iG}.Name, tru1.speed));
    fbd_save(f1, outDir, 'sensor_dynamics_1_structure');

    chk(end+1) = fbd_check('structure: wavelength recovered', an.lamPeak*1e3, lamP(iG)*1e3, ...
        0.15*lamP(iG)*1e3, 'mm');
    chk(end+1) = fbd_check('structure: frequency recovered', an.fPeak, f0, 1.0, 'Hz');
    chk(end+1) = fbd_check('structure: speed recovered', an.speed, tru1.speed, ...
        0.2*tru1.speed, 'm/s');

    % ---------- Fig 2: the same pattern through four instruments ----------
    f2 = fbd_fig(doExport, [60 60 1400 640]);
    for i = 1:nA
        [Ui, trui] = i_travelling(arrays{i}, pl{i}, lamP(i), f0, tv);
        subplot(2,nA,i);
        i_scatterfield(arrays{i}, Ui(:,1), max(abs(Ui(:))));
        title(sprintf('%s\n\\lambda = %.3g mm', arrays{i}.Name, lamP(i)*1e3), ...
            'Interpreter','tex', 'FontWeight','normal');
        subplot(2,nA,nA+i);
        i_carpet(Ui, pl{i}, arrays{i}, trui.khat, tv);
        if i == 1
            ylabel('sensor, ordered along k');
        end
        xlabel('time (s)');
    end
    sgtitle(['sensor\_dynamics -- Fig 2: one kind of wave, four instruments. A probe sees ' ...
             'a line; a grid sees a front; the helmet sees it curve with the surface']);
    fbd_save(f2, outDir, 'sensor_dynamics_2_instruments');

    % ---------- Fig 3: four dynamical families ----------
    fam = i_families(arrays{iG}, pl{iG}, G{iG}, cal{iG}, lamP(iG), f0, tv);
    f3 = fbd_fig(doExport, [70 70 1400 780]);
    % ⚠ NOT FRAME 1. A diffusing field starts as a delta and has spread nowhere yet, so the
    % first frame shows one lit sensor and says nothing. An eighth of the record in, every
    % family has something to show.
    jSnap = max(2, round(numel(tv)/8));
    for j = 1:numel(fam)
        subplot(3,numel(fam),j);
        i_scatterfield(arrays{iG}, fam(j).U(:,jSnap), max(abs(fam(j).U(:,jSnap))));
        title(fam(j).name, 'FontWeight','normal');
        subplot(3,numel(fam),numel(fam)+j);
        i_carpet(fam(j).U, pl{iG}, arrays{iG}, fam(j).khat, tv);
        if j==1, ylabel('sensor'); end
        subplot(3,numel(fam),2*numel(fam)+j);
        i_jointimage(fam(j).an);
        if j==1, ylabel('f (Hz)'); end
        xlabel('k (rad/m)');
    end
    sgtitle(['sensor\_dynamics -- Fig 3: four generators. Row 1 a snapshot, row 2 the ' ...
             'carpet, row 3 where each puts its energy in (k, f)']);
    fbd_save(f3, outDir, 'sensor_dynamics_3_families');

    % ---------- Fig 4: travelling vs standing -- the claim ----------
    [Utr, ~]  = i_travelling(arrays{iG}, pl{iG}, lamP(iG), f0, tv);
    Ust       = i_standing(arrays{iG}, pl{iG}, lamP(iG), f0, tv);
    ztr = i_atfreq(Utr, fs, f0);   zst = i_atfreq(Ust, fs, f0);
    sTr = i_travelindex(ztr);      sSt = i_travelindex(zst);
    iSensor = round(arrays{iG}.nCh/2);

    f4 = fbd_fig(doExport, [80 80 1300 720]);
    subplot(2,2,1);
    % ⚠ NORMALISED, AND THAT IS THE POINT. A single channel gives a sinusoid at f0 either
    % way. Its amplitude and phase depend only on WHERE that sensor happens to sit -- under a
    % standing wave's antinode or near its node -- and one channel has no reference against
    % which to judge either. So neither carries information about travelling vs standing.
    nrm = @(u) u / max(abs(u));
    plot(tv, nrm(Utr(iSensor,:)), '-', 'LineWidth',1.5); hold on;
    plot(tv, nrm(Ust(iSensor,:)), '--', 'LineWidth',1.5); grid on; xlim([0 0.4]);
    xlabel('time (s)'); ylabel('amplitude (normalised)');
    legend({'travelling','standing'}, 'Location','southeast');
    title(sprintf('one sensor (#%d): a %.0f Hz sinusoid, either way', iSensor, f0));
    subplot(2,2,2);
    fax = (0:nT-1)*(fs/nT);  keep = fax <= 40;
    ptr = abs(fft(Utr(iSensor,:) - mean(Utr(iSensor,:)))).^2;
    pst = abs(fft(Ust(iSensor,:) - mean(Ust(iSensor,:)))).^2;
    semilogy(fax(keep), ptr(keep), '-',  'LineWidth',1.3); hold on;
    semilogy(fax(keep), pst(keep), '--', 'LineWidth',1.3); grid on;
    xlabel('f (Hz)'); ylabel('power'); title('and the same single-channel spectrum');
    subplot(2,2,3);
    i_phasescatter(ztr);
    title(sprintf('travelling: phase advances with position  (index %.2f)', sTr));
    subplot(2,2,4);
    i_phasescatter(zst);
    title(sprintf('standing: every sensor on ONE line  (index %.2f)', sSt));
    sgtitle(['sensor\_dynamics -- Fig 4: indistinguishable per channel, opposite across the ' ...
             'array. This is what a graph-spectral method sees and a per-channel one cannot']);
    fbd_save(f4, outDir, 'sensor_dynamics_4_standing');

    chk(end+1) = fbd_check('travelling index (0.5 = travelling)', sTr, 0.5, 0.12, '');
    chk(end+1) = fbd_check('standing index (0 = standing)',       sSt, 0.0, 0.05, '');

    % ---------- Fig 5: a wave train of one speed -> a ridge -> the speed ----------
    aBig = rheome.sensors.grid('Size',[24 24], 'Pitch', 10e-3);
    Gb = rheome.sensors.graph(aBig);  cb = rheome.sensors.calibrate(Gb);  plb = sen_plane(aBig);
    cTrue = 0.35;                                            % m/s across the array
    [Utrain, kList, fList] = i_wavetrain(aBig, plb, cb, cTrue, tv);
    ab = i_analyse(Utrain, Gb, cb, fs);
    fit = i_ridgefit(ab);

    f5 = fbd_fig(doExport, [90 90 1300 480]);
    subplot(1,3,1);
    i_carpet(Utrain, plb, aBig, plb.e1.', tv);
    title('a train of components sharing one speed'); ylabel('sensor along k');
    subplot(1,3,2);
    i_jointimage(ab); hold on;
    kk = linspace(0, max(ab.kAxis), 50);
    plot(kk, fit.slope*kk/(2*pi), 'w--', 'LineWidth', 1.6);
    title('the ridge: \omega = c k');
    subplot(1,3,3);
    plot(kList, fList, 'o', 'MarkerSize', 7, 'LineWidth', 1.2); hold on;
    plot(kk, cTrue*kk/(2*pi), 'k-');
    plot(kk, fit.slope*kk/(2*pi), '--', 'LineWidth', 1.4); grid on;
    xlabel('k (rad/m)'); ylabel('f (Hz)');
    legend({'planted components','planted c','fitted c'}, 'Location','northwest');
    title(sprintf('c planted %.3f, fitted %.3f m/s (R^2 = %.3f)', cTrue, fit.slope, fit.r2));
    sgtitle(['sensor\_dynamics -- Fig 5: the dispersion relation is the measurement. ' ...
             'A straight ridge through the origin IS a non-dispersive travelling wave']);
    fbd_save(f5, outDir, 'sensor_dynamics_5_dispersion');

    chk(end+1) = fbd_check('wave-train speed recovered', fit.slope, cTrue, 0.15*cTrue, 'm/s');
    chk(end+1) = fbd_check('ridge fit quality', fit.r2, 1.0, 0.1, 'R^2');

    [ok, ~] = fbd_report(chk);

    if doExport
        i_writejson(outDir, arrays, pl, lamP, f0, tv, fs, iG);
    end

    out = struct('ok', ok, 'checks', chk, 'structure', an, 'families', fam, ...
                 'standingIndex', [sTr sSt], 'speedFit', fit);
end

% ===================== helpers =====================

function arr = i_meg()
    arr = [];
    try, arr = rheome.sensors.meg('subject01'); catch
        fprintf('  (subject01 not cached -- the MEG array is skipped)\n');
    end
end


function lam = i_pickwavelength(arr, c)
% ⚠ A WAVE THE ARRAY CAN SEE. Geometric centre of the usable band, so the demo measures
% dynamics rather than re-measuring the window rheome.demos.sensor_limits already reports.
    lo = c.lambdaUsable;
    if isnan(lo), lo = 4*arr.Pitch; end
    lam = sqrt(lo * arr.Aperture);
end


function [U, tru] = i_travelling(arr, pl, lam, f0, tv)
    k = 2*pi/lam;  w = 2*pi*f0;  khat = sen_dir(pl);
    x = pl.Pc * khat(:);
    U = rheome.sensors.sample(arr, @(P, t) cos(k*x - w*t), tv);
    tru = struct('k', k, 'khat', khat, 'f0', f0, 'lam', lam, 'speed', lam*f0);
end

function U = i_standing(arr, pl, lam, f0, tv)
    k = 2*pi/lam;  w = 2*pi*f0;  kh = sen_dir(pl);  x = pl.Pc * kh(:);
    U = rheome.sensors.sample(arr, @(P, t) cos(k*x) .* cos(w*t), tv);
end


function fam = i_families(arr, pl, G, cal, lam, f0, tv)
    k = 2*pi/lam;  w = 2*pi*f0;
    phi = atan2(pl.x2, pl.x1);  r = hypot(pl.x1, pl.x2);
    M = rheome.sensors.modes(G);
    d = zeros(arr.nCh,1);  [~,ic] = min(r);  d(ic) = 1;      % a delta at the centre
    Dif = M.Phi * (exp(-(M.Lambda/max(M.Lambda))*8*(tv)/max(tv)) .* (M.Phi.'*d));
    fam = struct('name', {}, 'U', {}, 'khat', {}, 'an', {});
    fam(1) = struct('name','travelling', 'U', sen_pattern(arr,pl,'travelling',lam,f0,tv), 'khat', pl.e1.', 'an', []);
    fam(2) = struct('name','rotating',   'U', sen_pattern(arr,pl,'rotating',  lam,f0,tv), 'khat', pl.e1.', 'an', []);
    fam(3) = struct('name','standing',   'U', sen_pattern(arr,pl,'standing',  lam,f0,tv), 'khat', pl.e1.', 'an', []);
    fam(4) = struct('name','diffusing',  'U', Dif, 'khat', pl.e1.', 'an', []);
    for j = 1:numel(fam)
        % Pass the calibration: without it the axis is sqrt(lambda), not rad/m, and
        % labelling a dimensionless axis in physical units is the error this whole
        % package exists to avoid.
        fam(j).an = i_analyse(fam(j).U, G, cal, 200);
    end
end

function an = i_analyse(U, G, cal, fs)
% Project onto graph modes, transform in time -> energy on the (k, f) plane.
    M  = rheome.sensors.modes(G);
    C  = M.Phi.' * (U - mean(U,2));
    nT = size(U,2);
    Cw = fft(C, [], 2) / sqrt(nT);
    nf = floor(nT/2) + 1;
    P  = abs(Cw(:,1:nf)).^2;
    fA = (0:nf-1) * (fs/nT);
    if isempty(cal), alpha = 1; else, alpha = cal.alpha; end
    kA = sqrt(max(M.Lambda,0) / alpha);                 % rad/m once alpha is known
    tot = sum(P(:)) + eps;
    [~, imax] = max(P(:));  [ik, jf] = ind2sub(size(P), imax);
    an = struct('P', P, 'kAxis', kA, 'fAxis', fA, ...
                'kPeak', kA(ik), 'fPeak', fA(jf), 'lamPeak', 2*pi/max(kA(ik),eps), ...
                'speed', 2*pi*fA(jf)/max(kA(ik),eps), ...
                'kBar', sum(kA.*sum(P,2))/tot, 'alpha', alpha);
end

function z = i_atfreq(U, fs, f0)
% The complex spatial pattern at one temporal frequency -- the analytic field at f0.
    nT = size(U,2);
    Uw = fft(U - mean(U,2), [], 2);
    fA = (0:nT-1)*(fs/nT);
    [~, j] = min(abs(fA(1:floor(nT/2)) - f0));
    z = Uw(:, j);
end

function i_phasescatter(z)
% The complex spatial pattern, with its principal axis drawn so collinearity is visible
% rather than inferred. axis equal, because a line must look like a line.
    A = [real(z(:)), imag(z(:))];  A = A - mean(A,1);
    [V, D] = eig(cov(A));  [~,ix] = max(diag(D));  v = V(:,ix);
    r = max(vecnorm(A,2,2));
    plot(A(:,1), A(:,2), '.', 'MarkerSize', 12); hold on;
    plot(r*[-v(1) v(1)], r*[-v(2) v(2)], '-', 'Color', [.55 .55 .6], 'LineWidth', 1.1);
    axis equal; grid on; xlabel('Re z'); ylabel('Im z');
end

function s = i_travelindex(z)
% ⭐ THE DISCRIMINATOR. Standing means one common phase, so every sensor lies on ONE line
% through the origin of the complex plane and the smaller eigenvalue of cov([Re z, Im z])
% vanishes. Travelling spreads the phase with position, filling the plane. The index is the
% eigenvalue split: 0 = standing, 0.5 = maximally travelling. It is amplitude-invariant.
    A = [real(z(:)), imag(z(:))];
    A = A - mean(A,1);
    e = sort(eig(cov(A)), 'descend');
    s = e(2) / max(sum(e), eps);
end

function [U, kList, fList] = i_wavetrain(arr, pl, c, cTrue, tv)
% Several plane waves that share ONE speed: f_j = c*k_j/(2*pi). Their joint spectrum is a
% straight ridge through the origin, and its slope is the speed.
    lamU = c.lambdaUsable;  if isnan(lamU), lamU = 4*arr.Pitch; end
    lams  = lamU * [1.15 1.5 2.0 2.7 3.6];
    lams  = lams(lams <= arr.Aperture);
    kList = 2*pi ./ lams;
    fList = cTrue * kList / (2*pi);
    U = zeros(arr.nCh, numel(tv));
    for j = 1:numel(kList)
        U = U + rheome.sensors.sample(arr, ...
            @(P,t) cos(kList(j)*pl.x1 - 2*pi*fList(j)*t + 0.7*j), tv);
    end
end

function fit = i_ridgefit(an)
% Energy-weighted fit of omega against k over the (k, f) energy -- the ridge's slope IS the
% speed. Weighted because bins carrying no signal outnumber those that do.
    % ⚠ THE MASK MUST EXCLUDE EMPTY BINS. A frequency carrying no signal still has a kbar --
    % the centroid of the whole k axis -- and those bins outnumber the occupied ones, so a
    % permissive threshold drags the fit toward a horizontal line. 2% of peak isolates the
    % components that are actually there.
    P = an.P;  w = sum(P, 1);  ok = w > 0.02*max(w);
    kbar = (an.kAxis(:).' * P(:,ok)) ./ w(ok);
    om   = 2*pi*an.fAxis(ok);
    wt   = (w(ok)/max(w(ok))).';
    A    = kbar(:);
    slope = (A.*wt).' * om(:) / max((A.*wt).' * A, eps);
    res  = om(:) - slope*A;
    ss   = sum(wt.*(om(:)-sum(wt.*om(:))/sum(wt)).^2);
    fit  = struct('slope', slope, 'r2', 1 - sum(wt.*res.^2)/max(ss,eps), ...
                  'kbar', kbar, 'f', an.fAxis(ok));
end

function i_scatterfield(arr, u, cl)
    P = arr.Pos * 1e3;
    if arr.Dim == 1
        scatter(P(:,3), zeros(size(P,1),1), 26, u, 'filled');
        ylim([-1 1]); xlabel('mm');
    else
        pl = sen_plane(arr);
        scatter(pl.x1*1e3, pl.x2*1e3, 26, u, 'filled');
    end
    axis equal tight off; caxis([-cl cl]); colormap(gca, i_diverging());
end

function i_traces(U, pl, arr, tv, n)
    [~, ord] = sort(pl.x1);
    pick = ord(round(linspace(1, numel(ord), n)));
    off = 2.2 * max(abs(U(:)));
    hold on;
    for j = 1:n
        plot(tv, U(pick(j),:) + (n-j)*off, 'LineWidth', 1.2);
    end
    grid on; xlim([0 max(tv)]); xlabel('time (s)');
    set(gca, 'YTick', (0:n-1)*off, 'YTickLabel', ...
        arrayfun(@(v) sprintf('%.3g mm', v*1e3), pl.x1(pick(end:-1:1)), 'UniformOutput', false));
end

function i_carpet(U, pl, arr, khat, tv)
    proj = pl.Pc * khat(:);
    [~, ord] = sort(proj);
    imagesc(tv, 1:size(U,1), U(ord,:));
    set(gca, 'YDir', 'normal'); caxis(max(abs(U(:)))*[-1 1]);
    colormap(gca, i_diverging()); xlabel('time (s)');
end

function i_jointimage(an)
% ⚠ CLIP TO THE TOP TWO DECADES. A full log scale renders the noise floor at nearly the same
% colour as the peak, so a spectrum that IS concentrated looks like a uniform wash.
    % ⚠ pcolor, NOT imagesc. The k axis is sqrt(lambda/alpha) and is NOT uniformly spaced;
    % imagesc ignores a non-uniform coordinate vector and lays cells out by INDEX between the
    % first and last value, so every feature is drawn at the wrong wavenumber. pcolor honours
    % the coordinates.
    L = log10(an.P.' + eps);
    top = max(L(:));
    pcolor(an.kAxis(:).', an.fAxis(:), L); shading flat;
    axis tight; ylim([0 min(40, max(an.fAxis))]);
    caxis([top-2 top]); colormap(gca, parula);
    xlabel('k (rad/m)'); ylabel('f (Hz)');
end

function m = i_diverging()
% Cool-to-warm through a light neutral: signed fields read as sign at a glance.
    n = 128; t = linspace(0,1,n).';
    lo = [0.13 0.29 0.55] + t.*([0.97 0.97 0.97] - [0.13 0.29 0.55]);
    hi = [0.97 0.97 0.97] + t.*([0.70 0.15 0.16] - [0.97 0.97 0.97]);
    m  = [lo; hi];
end

function i_writejson(outDir, arrays, pl, lamP, f0, tv, fs, iG)
% A compact space-time field per array, for an interactive rendering outside MATLAB.
    if ~exist(outDir,'dir'), mkdir(outDir); end
    % ⚠ FULL RATE, SHORTER RECORD. Decimating gave 10 frames per cycle at 10 Hz, which is
    % too coarse to watch a single cycle evolve. Keeping every frame over half the record
    % gives 20 frames per cycle for the same payload.
    keep = 1:min(128, numel(tv));
    S = struct('fs', fs, 'f0', f0, 't', tv(keep), 'arrays', []);
    % ⚠ THE INITIALISER GOVERNS. Assigning a struct with extra fields into a struct array
    % built with a shorter field list drops them silently -- no error, just missing data.
    A = struct('name', {}, 'kind', {}, 'x', {}, 'y', {}, 'dim', {}, 'lambda_mm', {}, ...
               'proj', {}, 'labels', {}, 'patterns', {}, 'fields', {});
    for i = 1:numel(arrays)
        names = {'travelling','standing','rotating','spiral','target','colliding'};
        F = struct();  def = {};
        for q = 1:numel(names)
            [Uq, okq] = sen_pattern(arrays{i}, pl{i}, names{q}, lamP(i), f0, tv);
            if ~okq, continue; end
            sq = max(abs(Uq(:)));  if sq <= 0, sq = 1; end
            F.(names{q}) = round(Uq(:,keep)/sq, 2);
            def{end+1} = names{q}; %#ok<AGROW>
        end
        % Carry the planting direction and the channel names: a viewer needs to order the
        % traces ALONG propagation to see the lag, and to name the channel it is looking at.
        kh = sen_dir(pl{i});
        A(end+1) = struct('name', arrays{i}.Name, 'kind', arrays{i}.Kind, ...
            'x', round(pl{i}.x1(:).'*1e3, 3), 'y', round(pl{i}.x2(:).'*1e3, 3), ...
            'dim', arrays{i}.Dim, 'lambda_mm', round(lamP(i)*1e3, 3), ...
            'proj', round((pl{i}.Pc * kh(:)).'*1e3, 3), ...
            'labels', {arrays{i}.Labels}, ...
            'patterns', {def}, 'fields', F); %#ok<AGROW>
    end
    S.arrays = A;
    fid = fopen(fullfile(outDir, 'sensor_dynamics.json'), 'w');
    fwrite(fid, jsonencode(S));  fclose(fid);
    fprintf('  wrote %s\n', fullfile(outDir, 'sensor_dynamics.json'));
end

% Author: Diellor Basha, 2026
