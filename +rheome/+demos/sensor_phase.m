function out = sensor_phase(outDir)
% DEMOS.SENSOR_PHASE  Instantaneous phase across an array, and what its gradient measures.
%
%   rheome.demos.sensor_phase()
%   out = rheome.demos.sensor_phase('_figures')
%
%   1 maps      amplitude and phase side by side: travelling ramps, standing steps by pi
%   2 montage   unwrapped phase against position -- a straight line whose slope IS k
%   3 gradient  rheome.flow.phasegradient on the sensor graph: a wavevector per face
%   4 rotor     a rotating wave, and rheome.detect.phasesingularity finding its core and chirality
%   5 lag       phase lag against sensor separation -- and why PLV alone cannot see travel
%
% ⭐ NOTHING HERE IS NEW MACHINERY. rheome.flow.phasegradient and rheome.detect.phasesingularity were
% written for triangulated surface meshes and are used UNCHANGED, because rheome.sensors.graph
% returns the same .Vertices/.Faces/.nV shape a mesh does. That was the point of
% building it that way.
%
% ⭐ THE ANALYTIC SIGNAL COMES FROM THE FFT, NOT FROM hilbert(). Zeroing the negative
% frequencies and doubling the positive ones IS the Hilbert transform, it needs no toolbox,
% and it is the same route rheome.flow.joint takes.
%
% ⚠ PHASE IS THE MEASUREMENT, AMPLITUDE IS NOT. Multiply the field by any positive envelope
% -- a gain change, a source waxing and waning -- and every phase quantity here is untouched
% while every amplitude one moves. That is why a phase gradient is a more stable readout of
% propagation than a peak-tracking one.
%
% ⚠ SENSOR SPACE ONLY. k is rad per metre ACROSS THE ARRAY; speed is metres per second
% across it. Nothing here is a claim about what produced the signal.
%
% INPUTS:
%   outDir  '' for screen, or a folder for PNGs
% OUTPUT:
%   out  .ok .checks .kFromPhase .kFromGradient .rotor .lagSlope
%
% See also: rheome.demos.sensor_dynamics, rheome.flow.phasegradient, rheome.detect.phasesingularity
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});
    fs = 200;  nT = 256;  tv = (0:nT-1)/fs;  f0 = 10;

    arr = rheome.sensors.grid('ecog');
    G   = rheome.sensors.graph(arr);
    cal = rheome.sensors.calibrate(G);
    pl  = i_plane(arr);
    lam = sqrt(cal.lambdaUsable * arr.Aperture);
    kTrue = 2*pi/lam;  khat = i_dir(pl);
    fprintf('\n=== rheome.demos.sensor_phase : %s, %.1f mm wave at %g Hz ===\n', arr.Name, lam*1e3, f0);

    Utr = i_travelling(arr, pl, lam, f0, tv);
    Ust = i_standing(arr, pl, lam, f0, tv);
    ztr = i_analytic(Utr);   zst = i_analytic(Ust);
    jm  = round(nT/2);                                   % a frame away from the edges

    % ---------- Fig 1: amplitude and phase, travelling vs standing ----------
    f1 = fbd_fig(doExport, [50 50 1300 660]);
    pans = {abs(ztr(:,jm)), angle(ztr(:,jm)), abs(zst(:,jm)), angle(zst(:,jm))};
    ttl  = {'travelling: |z|', 'travelling: phase', 'standing: |z|', 'standing: phase'};
    for j = 1:4
        subplot(2,2,j);
        isPhase = mod(j,2)==0;
        i_scatter(arr, pl, pans{j}, isPhase);
        title(ttl{j}, 'FontWeight','normal'); colorbar;
    end
    sgtitle(['sensor\_phase -- Fig 1: a travelling wave RAMPS its phase across the array; ' ...
             'a standing wave holds it constant and flips by \pi at each node']);
    fbd_save(f1, outDir, 'sensor_phase_1_maps');

    % ---------- Fig 2: phase along a montage -> the slope IS k ----------
    mont = i_montage(pl, khat, arr);
    xm   = pl.Pc(mont,:) * khat(:);
    phTr = unwrap(angle(ztr(mont, jm)));
    phSt = unwrap(angle(zst(mont, jm)));
    pf   = polyfit(xm, phTr, 1);   kPhase = abs(pf(1));
    pfS  = polyfit(xm, phSt, 1);

    f2 = fbd_fig(doExport, [60 60 1250 460]);
    subplot(1,2,1);
    plot(xm*1e3, phTr, 'o-', 'LineWidth',1.3, 'MarkerSize',5); hold on;
    plot(xm*1e3, polyval(pf, xm), 'k--'); grid on;
    xlabel('position along k (mm)'); ylabel('unwrapped phase (rad)');
    title(sprintf('travelling: slope %.0f rad/m -> \\lambda = %.1f mm (planted %.1f)', ...
        kPhase, 2*pi/kPhase*1e3, lam*1e3));
    subplot(1,2,2);
    plot(xm*1e3, phSt, 'o-', 'LineWidth',1.3, 'MarkerSize',5); grid on;
    xlabel('position along k (mm)'); ylabel('unwrapped phase (rad)');
    title(sprintf('standing: slope %.0f rad/m -- flat, with \\pi steps', abs(pfS(1))));
    sgtitle(['sensor\_phase -- Fig 2: phase against position. A constant slope is a wave ' ...
             'going somewhere; a flat line with jumps is one going nowhere']);
    fbd_save(f2, outDir, 'sensor_phase_2_montage');

    chk(end+1) = fbd_check('wavelength from phase slope', 2*pi/kPhase*1e3, lam*1e3, ...
        0.12*lam*1e3, 'mm');
    chk(end+1) = fbd_check('standing slope is ~zero', abs(pfS(1)), 0, 0.15*kTrue, 'rad/m');

    % ---------- Fig 3: the phase gradient, per face ----------
    S  = struct('Vertices', G.Vertices, 'Faces', G.Faces, 'nV', G.nV);
    pg = rheome.flow.phasegradient(ztr(:, jm-2:jm+2), S, 'Rate', fs);
    kmag = median(pg.kmag(:,3), 'omitnan');
    spd  = median(pg.speed(:,3), 'omitnan');

    f3 = fbd_fig(doExport, [70 70 1300 460]);
    subplot(1,3,1);
    i_quiver(S, pl, pg, 3);
    title('wavevector per face');
    subplot(1,3,2);
    histogram(pg.kmag(:,3), 24); grid on; hold on;
    xline(kTrue, 'r-', 'LineWidth', 1.5);
    xlabel('|k| (rad/m)'); ylabel('faces');
    title(sprintf('|k| = %.0f, planted %.0f rad/m', kmag, kTrue));
    subplot(1,3,3);
    histogram(pg.speed(:,3), 24); grid on; hold on;
    xline(lam*f0, 'r-', 'LineWidth', 1.5);
    xlabel('speed (m/s)'); ylabel('faces');
    title(sprintf('speed %.2f, planted %.2f m/s', spd, lam*f0));
    sgtitle(sprintf(['sensor\\_phase -- Fig 3: rheome.flow.phasegradient, unchanged from the mesh ' ...
        'code. Aliased faces: %.0f%%'], 100*pg.aliased(3)));
    fbd_save(f3, outDir, 'sensor_phase_3_gradient');

    chk(end+1) = fbd_check('|k| from phase gradient', kmag, kTrue, 0.12*kTrue, 'rad/m');
    chk(end+1) = fbd_check('speed from phase gradient', spd, lam*f0, 0.15*lam*f0, 'm/s');
    chk(end+1) = fbd_check('aliased faces', pg.aliased(3), 0, 0.05, 'frac');

    % ---------- Fig 4: a rotating wave and its singularity ----------
    phi = atan2(pl.x2, pl.x1);
    Urot = rheome.sensors.sample(arr, @(P,t) cos(phi - 2*pi*f0*t), tv);
    zrot = i_analytic(Urot);
    ps   = rheome.detect.phasesingularity(zrot(:, jm), S);
    nCore = numel(ps.charge);
    coreErr = NaN;
    if nCore >= 1
        cen = mean(arr.Pos, 1);
        d = vecnorm(ps.pos - cen, 2, 2);
        [coreErr, ib] = min(d);
    end

    f4 = fbd_fig(doExport, [80 80 1250 460]);
    subplot(1,2,1);
    i_scatter(arr, pl, angle(zrot(:,jm)), true); colorbar; hold on;
    if nCore >= 1
        pc = (ps.pos(ib,:) - mean(arr.Pos,1)) * [pl.e1 pl.e2];
        plot(pc(1)*1e3, pc(2)*1e3, 'ko', 'MarkerSize', 13, 'LineWidth', 2);
    end
    title('phase of a rotating wave, with the detected core');
    subplot(1,2,2);
    i_faceimage(S, pl, ps.perFace);
    title(sprintf('winding per face: %d core(s), charge %s', nCore, ...
        mat2str(unique(ps.charge(:)).')));
    sgtitle(['sensor\_phase -- Fig 4: the phase winds once around a point. An integer ' ...
             'cannot drift, which is why a winding number beats a threshold on curl']);
    fbd_save(f4, outDir, 'sensor_phase_4_rotor');

    chk(end+1) = fbd_check('rotor: one core found', nCore, 1, 0.5, '');
    chk(end+1) = fbd_check('rotor: core near array centre', coreErr*1e3, 0, ...
        0.25*arr.Aperture*1e3, 'mm');

    % ---------- Fig 5: phase lag against separation, and the PLV trap ----------
    ref  = 1;
    dsep = pl.Pc * khat(:);  dsep = dsep - dsep(ref);
    dphi = angle(ztr(:,jm) .* conj(ztr(ref,jm)));
    dphiU = unwrap_by(dsep, dphi);
    pl2  = polyfit(dsep, dphiU, 1);   kLag = abs(pl2(1));
    plv  = abs(mean(exp(1i*angle(ztr .* conj(ztr(ref,:)))), 2));

    f5 = fbd_fig(doExport, [90 90 1250 460]);
    subplot(1,2,1);
    plot(dsep*1e3, dphiU, '.', 'MarkerSize', 12); hold on;
    plot(dsep*1e3, polyval(pl2, dsep), 'k--'); grid on;
    xlabel('separation along k (mm)'); ylabel('phase lag vs reference (rad)');
    title(sprintf('lag grows with distance: %.0f rad/m -> %.1f mm', kLag, 2*pi/kLag*1e3));
    subplot(1,2,2);
    plot(abs(dsep)*1e3, plv, '.', 'MarkerSize', 12); grid on; ylim([0 1.05]);
    xlabel('separation (mm)'); ylabel('PLV vs reference');
    title('PLV is 1 everywhere -- it cannot see the travel');
    sgtitle(['sensor\_phase -- Fig 5: two channels perfectly phase-locked can still be ' ...
             'a travelling wave. The information is in the LAG, not the locking']);
    fbd_save(f5, outDir, 'sensor_phase_5_lag');

    chk(end+1) = fbd_check('wavelength from pairwise lag', 2*pi/kLag*1e3, lam*1e3, ...
        0.12*lam*1e3, 'mm');
    chk(end+1) = fbd_check('PLV is ~1 for every pair', min(plv), 1, 0.05, '');

    [ok, ~] = fbd_report(chk);
    out = struct('ok', ok, 'checks', chk, 'kFromPhase', kPhase, 'kFromGradient', kmag, ...
                 'rotor', ps, 'lagSlope', kLag);
end

% ===================== helpers =====================

function z = i_analytic(U)
% ⭐ THE ANALYTIC SIGNAL WITHOUT A TOOLBOX. Zero the negative frequencies, double the
% positive ones: that IS the Hilbert transform, and it is what rheome.flow.joint does.
    U = U - mean(U, 2);
    nT = size(U,2);
    F  = fft(U, [], 2);
    h  = zeros(1, nT);
    h(1) = 1;
    if mod(nT,2)==0, h(nT/2+1) = 1; h(2:nT/2) = 2; else, h(2:(nT+1)/2) = 2; end
    z  = ifft(F .* h, [], 2);
end

function pl = i_plane(arr)
    Pc = arr.Pos - mean(arr.Pos, 1);
    [~,~,V] = svd(Pc, 'econ');
    if arr.Dim < 2
        e1 = V(:,1); e2 = V(:,2);
    else
        [~, drop] = max(abs(V(:,3)));  keepAx = setdiff(1:3, drop);
        E = eye(3);  e1 = E(:,keepAx(1));  e2 = E(:,keepAx(2));
    end
    pl = struct('e1', e1, 'e2', e2, 'x1', Pc*e1, 'x2', Pc*e2, 'Pc', Pc);
end

function d = i_dir(pl)
    d = (cosd(30)*pl.e1 + sind(30)*pl.e2).';  d = d / norm(d);
end

function U = i_travelling(arr, pl, lam, f0, tv)
    kh = i_dir(pl);  x = pl.Pc * kh(:);
    U = rheome.sensors.sample(arr, @(P,t) cos(2*pi/lam*x - 2*pi*f0*t), tv);
end

function U = i_standing(arr, pl, lam, f0, tv)
    kh = i_dir(pl);  x = pl.Pc * kh(:);
    U = rheome.sensors.sample(arr, @(P,t) cos(2*pi/lam*x) .* cos(2*pi*f0*t), tv);
end

function idx = i_montage(pl, khat, arr)
% A LINE of sensors through the array, ordered along it -- neighbouring traces are
% neighbouring sensors, which is what makes a phase-versus-position plot mean anything.
    perp = [-khat(2) khat(1) 0];
    t = pl.Pc * khat(:);  u = pl.Pc * perp(:);
    band = 0.08 * (max(u) - min(u) + eps);
    idx = find(abs(u - median(u)) <= band);
    if numel(idx) < 4, [~, o] = sort(abs(u - median(u))); idx = o(1:min(8,numel(o))); end
    [~, o2] = sort(t(idx));  idx = idx(o2);
end

function y = unwrap_by(x, phi)
% Unwrap a phase sampled at scattered positions: sort by position, unwrap, put back.
    [~, o] = sort(x);
    y = zeros(size(phi));
    y(o) = unwrap(phi(o));
end

function i_scatter(arr, pl, v, isPhase)
    scatter(pl.x1*1e3, pl.x2*1e3, 46, v, 'filled');
    axis equal tight; xlabel('mm');
    if isPhase
        colormap(gca, hsv(256)); caxis([-pi pi]);
    else
        colormap(gca, parula); caxis([0 max(v)]);
    end
end

function i_quiver(S, pl, pg, j)
    C = (S.Vertices(S.Faces(:,1),:) + S.Vertices(S.Faces(:,2),:) + S.Vertices(S.Faces(:,3),:))/3;
    Cc = C - mean(S.Vertices,1);
    kx = pg.k(:,:,j) * pl.e1;  ky = pg.k(:,:,j) * pl.e2;
    quiver(Cc*pl.e1*1e3, Cc*pl.e2*1e3, kx, ky, 1.1, 'LineWidth', 1);
    axis equal tight; grid on; xlabel('mm'); ylabel('mm');
end

function i_faceimage(S, pl, perFace)
    C = (S.Vertices(S.Faces(:,1),:) + S.Vertices(S.Faces(:,2),:) + S.Vertices(S.Faces(:,3),:))/3;
    Cc = C - mean(S.Vertices,1);
    scatter(Cc*pl.e1*1e3, Cc*pl.e2*1e3, 30, perFace, 'filled');
    axis equal tight; colorbar; caxis([-1 1]);
    colormap(gca, [0 0.28 0.62; 0.93 0.93 0.93; 0.70 0.15 0.16]);
    xlabel('mm'); ylabel('mm');
end

% Author: Diellor Basha, 2026
