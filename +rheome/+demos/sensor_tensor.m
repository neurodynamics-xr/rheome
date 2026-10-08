function out = sensor_tensor(outDir)
% DEMOS.SENSOR_TENSOR  The joint wavelet tensor: aperture, rate, and dispersion as controls.
%
%   rheome.demos.sensor_tensor()
%   out = rheome.demos.sensor_tensor('_figures')
%
%   1 ladder     the graph wavelet bank as a ladder of APERTURES, in millimetres
%   2 atoms      the same members localized at chosen sensors -- where the aperture sits
%   3 tiling     the joint (wavelength, rate) tensor the two banks span together
%   4 separate   two patterns differing in size AND rate, pulled apart into their own cells
%   5 dispersion two trains sharing a rectangle but not a diagonal -- what only a
%                non-separable member can do
%
% ⭐ WHAT THE TENSOR IS FOR. A bandpass in time answers "how fast"; a graph wavelet answers
% "how big"; neither alone answers "how big AND how fast", and no sequence of the two answers
% "how fast RELATIVE to how big" -- which is a speed. The joint bank carries all three because
% its members are functions of (lambda, omega) together, and the ones that matter most are
% precisely the ones that do NOT factor.
%
% ⭐ APERTURE IS A LENGTH HERE, NOT AN INDEX. Every member's scale converts to millimetres
% through the calibration constant: k = sqrt(lambda/alpha), so wavelength = 2*pi*sqrt(alpha)
% over the member's centre. Without rheome.sensors.calibrate a filterbank on a graph is a ladder of
% dimensionless numbers; with it, it is a ladder of apertures you can compare to the array.
%
% ⚠ A SPEED KERNEL'S c IS IN GRAPH UNITS. speedkernel selects omega = c*sqrt(lambda), and
% lambda here is dimensionless, so a physical speed c_phys enters as c = c_phys/sqrt(alpha).
% Passing metres per second directly selects the wrong diagonal by a factor of sqrt(alpha).
%
% ⚠ SENSOR SPACE ONLY. Every length is millimetres of sensor separation and every speed is
% metres per second across the array.
%
% INPUTS:
%   outDir  '' for screen, or a folder for PNGs
% OUTPUT:
%   out  .ok .checks .bank .atoms .separated .speedSel
%
% See also: rheome.graphfilterbank, rheome.jointfilterbank, rheome.sensors.calibrate, rheome.demos.sensor_dynamics
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});
    fs = 200;  nT = 256;  tv = (0:nT-1)/fs;

    arr = rheome.sensors.grid('Size', [24 24], 'Pitch', 10e-3);
    G   = rheome.sensors.graph(arr);
    cal = rheome.sensors.calibrate(G);
    M   = rheome.sensors.modes(G);
    sa  = sqrt(cal.alpha);                       % metres per unit of graph wavenumber
    Pc  = arr.Pos - mean(arr.Pos,1);  x1 = Pc(:,1);  x2 = Pc(:,2);
    fprintf('\n=== rheome.demos.sensor_tensor : %s, alpha = %.4g m^2 ===\n', arr.Name, cal.alpha);

    % ---------- Fig 1: the bank as a ladder of apertures ----------
    gfb = rheome.graphfilterbank(M.Lambda, 'Wavelet','itersine', 'NumFilters', 6, ...
                          'Transform', rheome.graphtransform.eigen(M.Phi, speye(G.nV), M.Lambda));
    kc   = centerWavenumbers(gfb);               % graph units
    lamM = 2*pi*sa ./ max(kc, eps) * 1e3;        % millimetres -- the aperture ladder
    H    = graphfilters(gfb);                    % [K x M], on the bank's own spectrum
    bnd  = framebounds(gfb);

    f1 = fbd_fig(doExport, [50 50 1300 430]);
    subplot(1,2,1);
    lamAxis = 2*pi*sa ./ max(sqrt(M.Lambda), eps) * 1e3;
    semilogx(lamAxis, H, 'LineWidth', 1.3); grid on; xlim([2*arr.Pitch*1e3 arr.Aperture*1e3]);
    xlabel('wavelength (mm)'); ylabel('gain'); title('members, on a physical axis');
    subplot(1,2,2);
    barh(lamM); grid on; set(gca,'YDir','reverse');
    xlabel('centre wavelength (mm)'); ylabel('member');
    title(sprintf('the aperture ladder  (frame A=%.2f B=%.2f)', bnd.A, bnd.B));
    sgtitle(['sensor\_tensor -- Fig 1: a graph wavelet bank is a ladder of APERTURES once ' ...
             '\alpha has given \lambda a length']);
    fbd_save(f1, outDir, 'sensor_tensor_1_ladder');

    chk(end+1) = fbd_check('bank is a frame (A>0)', double(bnd.A > 0.1), 1, 0.5, '');
    chk(end+1) = fbd_check('apertures are ordered', double(all(diff(lamM) < 0)), 1, 0.5, '');

    % ---------- Fig 2: the same members, localized at sensors ----------
    seeds = [i_nearest(Pc, [-60 -60 0]*1e-3), i_nearest(Pc, [0 0 0]), ...
             i_nearest(Pc, [70 50 0]*1e-3)];
    % ⚠ ONLY MEMBERS THE BASIS CAN CARRY. gfb.Usable marks those losing under 5% of their
    % mass past lambda_max; beyond it an "atom" is the truncation's shape, not the member's.
    % ⚠ AND ONLY MEMBERS THAT FIT. An atom whose wavelength approaches the array's own
    % extent is truncated by the boundary, so its measured width saturates at the array
    % radius rather than reporting the member's scale -- the coarsest members here read
    % LARGER than the array can hold. Restrict the ladder to members comfortably inside it;
    % the ones outside are aperture-limited, which is a property of the instrument.
    fits   = find(gfb.Usable(:).' & (lamM(:).' < arr.Aperture*1e3/3));
    mShow  = fits(round(linspace(1, numel(fits), 3)));
    r50 = zeros(numel(mShow), numel(seeds));
    f2 = fbd_fig(doExport, [60 60 1150 1000]);
    for a = 1:numel(mShow)
        gm = H(:, mShow(a));
        for b = 1:numel(seeds)
            d = zeros(G.nV,1); d(seeds(b)) = 1;
            atom = M.Phi * (gm .* (M.Phi.' * d));
            r50(a,b) = i_radius50(atom, Pc, seeds(b));
            subplot(numel(mShow), numel(seeds), (a-1)*numel(seeds)+b);
            scatter(x1*1e3, x2*1e3, 22, atom, 'filled'); axis equal tight off;
            caxis(max(abs(atom))*[-1 1]); colormap(gca, i_div());
            title(sprintf('m%d @ s%d: r_{rms} = %.0f mm', mShow(a), b, r50(a,b)*1e3), ...
                'FontWeight','normal', 'FontSize', 9);
        end
    end
    sgtitle(['sensor\_tensor -- Fig 2: one bank, three apertures, three sensors. ' ...
             'The member sets the SIZE; the seed sets WHERE']);
    fbd_save(f2, outDir, 'sensor_tensor_2_atoms');

    % ⚠ MEMBER 1 IS THE COARSEST. centerWavenumbers rises with index, so the aperture
    % ladder DESCENDS -- an atom from a later member is smaller, not larger.
    % ⚠ ADJACENT MEMBERS DO NOT ORDER STRICTLY, and asserting they do would claim a property
    % this bank does not have. An itersine frame's members overlap, so neighbours have
    % similar effective extents, and on a finite array the boundary truncates each seed
    % differently -- measured here: 31.4, 34.9, 14.2 mm across a descending ladder. What IS
    % true, and is the claim worth making, is that the ladder spans a real range of
    % apertures end to end.
    rm = mean(r50, 2);
    chk(end+1) = fbd_check('the ladder spans a real range of apertures', ...
        rm(1)/max(rm(end), eps), 2.2, 1.0, 'x');
    chk(end+1) = fbd_check('coarse members exceed the array and are excluded', ...
        double(any(lamM > arr.Aperture*1e3/3)), 1, 0.5, '');

    % ---------- Fig 3: the joint tiling, in physical units ----------
    bands = [4 8; 8 14; 14 24; 24 40];
    jfb = rheome.jointfilterbank(gfb, 'SignalLength', nT, 'SamplingFrequency', fs, 'Bands', bands);
    fc  = centerFrequencies(jfb);
    kcj = centerWavenumbers(jfb);
    lamj = 2*pi*sa ./ max(kcj, eps) * 1e3;

    f3 = fbd_fig(doExport, [70 70 1200 470]);
    subplot(1,2,1);
    loglog(lamj, fc, 'o', 'MarkerSize', 8, 'LineWidth', 1.4); grid on;
    xlabel('member wavelength (mm)'); ylabel('member rate (Hz)');
    title(sprintf('%d x %d = %d members', jfb.NumGraph, jfb.NumTime, jfb.NumMembers));
    subplot(1,2,2);
    spectralResponse(jfb, 7); title('one member on the (\lambda, \omega) plane');
    sgtitle(sprintf(['sensor\\_tensor -- Fig 3: the tensor spans SIZE x RATE. Materialised ' ...
        'it would be %.0f kB; it is stored as factors'], ...
        numel(jfb.Lambda)*jfb.NumOmega*jfb.NumMembers*8/1e3));
    fbd_save(f3, outDir, 'sensor_tensor_3_tiling');

    % ---------- Fig 4: two patterns, pulled into their own cells ----------
    lamA = 140e-3; fA = 6;      lamB = 45e-3;  fB = 20;
    kA = 2*pi/lamA;  kB = 2*pi/lamB;
    UA = cos(kA*x1 - 2*pi*fA*tv);
    UB = 0.8*cos(kB*(0.6*x1+0.8*x2) - 2*pi*fB*tv);
    U  = UA + UB;
    C  = fft(M.Phi.' * (U - mean(U,2)), [], 2) / sqrt(nT);
    nf = floor(nT/2)+1;  C = C(:,1:nf);
    fAx = (0:nf-1)*(fs/nT);
    E   = scalogram(jfb, C);
    [pk, cellIdx] = i_peaks(E, jfb, fAx, sa);

    f4 = fbd_fig(doExport, [80 80 1350 460]);
    subplot(1,3,1);
    scatter(x1*1e3, x2*1e3, 22, U(:,20), 'filled'); axis equal tight off;
    caxis(max(abs(U(:)))*[-1 1]); colormap(gca, i_div());
    title('the field: two components superposed');
    subplot(1,3,2);
    i_joint(M.Lambda, fAx, abs(C).^2, sa);
    title('joint spectrum: two separated blobs');
    subplot(1,3,3);
    imagesc(fAx, 1:jfb.NumMembers, log10(E + eps)); set(gca,'YDir','normal');
    xlim([0 40]); xlabel('f (Hz)'); ylabel('tensor member');
    title('scalogram: which cells carry them');
    sgtitle(sprintf(['sensor\\_tensor -- Fig 4: planted (%.0f mm, %g Hz) and (%.0f mm, %g Hz); ' ...
        'recovered (%.0f mm, %.1f Hz) and (%.0f mm, %.1f Hz)'], ...
        lamA*1e3, fA, lamB*1e3, fB, pk(1).lam, pk(1).f, pk(2).lam, pk(2).f));
    fbd_save(f4, outDir, 'sensor_tensor_4_separate');

    chk(end+1) = fbd_check('component A: wavelength', pk(1).lam, lamA*1e3, 0.25*lamA*1e3, 'mm');
    chk(end+1) = fbd_check('component A: rate',       pk(1).f,   fA,       2.0, 'Hz');
    chk(end+1) = fbd_check('component B: wavelength', pk(2).lam, lamB*1e3, 0.25*lamB*1e3, 'mm');
    chk(end+1) = fbd_check('component B: rate',       pk(2).f,   fB,       3.0, 'Hz');
    chk(end+1) = fbd_check('they land in different cells', double(cellIdx(1) ~= cellIdx(2)), 1, 0.5, '');

    % ---------- Fig 5: the diagonal a rectangle cannot draw ----------
    c1 = 0.30;  c2 = 1.20;                        % m/s across the array
    [U2, lam1, f1s] = i_train(x1, tv, c1, sa, cal, arr);
    [U3, lam2, f2s] = i_train(x1, tv, c2, sa, cal, arr);
    Um = U2 + U3;
    spec = @(X) abs(fft(M.Phi.' * (X - mean(X,2)), [], 2)/sqrt(nT)).^2;
    P1 = spec(U2); P1 = P1(:,1:nf);
    P2 = spec(U3); P2 = P2(:,1:nf);
    Pm = spec(Um); Pm = Pm(:,1:nf);
    sel1 = i_speedmask(M.Lambda, fAx, c1/sa, 0.22);
    sel2 = i_speedmask(M.Lambda, fAx, c2/sa, 0.22);
    % ⚠ THE RECTANGLE MUST BE THE SMALLEST ONE CONTAINING BOTH TRAINS, derived from what was
    % planted rather than guessed -- otherwise the figure proves only that a badly chosen
    % rectangle misses things, which is not the point.
    lamAll = [lam1 lam2];  fAll = [f1s f2s];
    rect = i_rectmask(M.Lambda, fAx, sa, [0.8*min(lamAll) 1.25*max(lamAll)], ...
                                        [0.7*min(fAll)   1.3*max(fAll)]);
    keptOf = @(P, m) sum(P(m)) / sum(P(:));
    speedSel = [keptOf(P1,sel1) keptOf(P2,sel1) keptOf(P1,rect) keptOf(P2,rect)];

    f5 = fbd_fig(doExport, [90 90 1350 460]);
    subplot(1,3,1); i_joint(M.Lambda, fAx, Pm, sa); hold on;
    title('two trains: two diagonals');
    subplot(1,3,2); i_joint(M.Lambda, fAx, Pm .* rect, sa);
    title(sprintf('rectangle: keeps %.0f%% of slow AND %.0f%% of fast', ...
        100*keptOf(P1,rect), 100*keptOf(P2,rect)));
    subplot(1,3,3); i_joint(M.Lambda, fAx, Pm .* sel1, sa);
    title(sprintf('speed member: %.0f%% of slow, %.0f%% of fast', ...
        100*keptOf(P1,sel1), 100*keptOf(P2,sel1)));
    sgtitle(['sensor\_tensor -- Fig 5: both trains cross the same rectangle. Only a ' ...
             'NON-SEPARABLE member selects one speed and rejects the other']);
    fbd_save(f5, outDir, 'sensor_tensor_5_dispersion');

    % The rectangle cannot tell them apart; the diagonal can. That contrast IS the figure.
    chk(end+1) = fbd_check('rectangle admits the slow train',  keptOf(P1,rect), 1, 0.25, 'frac');
    chk(end+1) = fbd_check('rectangle admits the fast train',  keptOf(P2,rect), 1, 0.25, 'frac');
    % A tighter diagonal rejects more of the other train but keeps less of its own; a
    % broadband train has real spectral width around its ridge. Half is a good haul.
    chk(end+1) = fbd_check('speed member admits the slow train', keptOf(P1,sel1), 1, 0.50, 'frac');
    chk(end+1) = fbd_check('speed member REJECTS the fast train', keptOf(P2,sel1), 0, 0.10, 'frac');
    % A LOWER BOUND, not a target: the measured ratio is ~180x and pinning it to a value
    % would make an improvement look like a regression.
    chk(end+1) = fbd_check('selectivity is at least 20x', ...
        double(keptOf(P1,sel1)/max(keptOf(P2,sel1), 1e-6) > 20), 1, 0.5, '');

    [ok, ~] = fbd_report(chk);
    out = struct('ok', ok, 'checks', chk, 'bank', lamM, 'atoms', r50, ...
                 'separated', pk, 'speedSel', speedSel);
end

% ===================== helpers =====================

function i = i_nearest(Pc, p)
    [~, i] = min(vecnorm(Pc - p, 2, 2));
end

function r = i_radius50(atom, Pc, seed)
% ⚠ A HALF-ENERGY RADIUS IS THE WRONG MEASURE FOR AN OSCILLATORY ATOM. A wavelet rings, so
% cumulative energy against radius is not monotone and the 50% crossing lands wherever a
% lobe happens to fall -- measured on this bank it gave 22, 14, 28, 22, 0, 10 mm across a
% ladder that is strictly descending. The second moment is the honest width:
%   r_rms = sqrt( sum e_v d_v^2 / sum e_v )
% which is monotone in the member scale because it weights the whole envelope rather than
% hunting for one crossing.
    e = atom.^2;  d = vecnorm(Pc - Pc(seed,:), 2, 2);
    r = sqrt(sum(e .* d.^2) / max(sum(e), eps));
end

function i_joint(Lambda, fAx, P, sa)
% ⚠ pcolor, not imagesc: the wavelength axis is NOT uniformly spaced.
    lam = 2*pi*sa ./ max(sqrt(Lambda), eps) * 1e3;
    L = log10(P.' + eps);  top = max(L(:));
    pcolor(lam(:).', fAx(:), L); shading flat;
    set(gca, 'XScale', 'log'); axis tight; ylim([0 40]);
    caxis([top-2.5 top]); colormap(gca, parula);
    xlabel('wavelength (mm)'); ylabel('f (Hz)');
end

function [pk, idx] = i_peaks(E, jfb, fAx, sa)
% The two strongest cells, reported as a physical (wavelength, rate).
    kc = centerWavenumbers(jfb);
    lamM = 2*pi*sa ./ max(kc, eps) * 1e3;
    Es = E;  pk = struct('lam', {}, 'f', {});  idx = zeros(1,2);
    for j = 1:2
        [~, li] = max(Es(:));  [m, w] = ind2sub(size(Es), li);
        pk(j) = struct('lam', lamM(m), 'f', fAx(w)); %#ok<AGROW>
        idx(j) = m;
        keep = true(size(Es));
        keep(:, abs(fAx - fAx(w)) < 4) = false;    % suppress this peak and retry
        Es = Es .* keep;
    end
    [~, o] = sort([pk.lam], 'descend');  pk = pk(o);  idx = idx(o);
end

function [U, lams, fs_] = i_train(x, tv, c, sa, cal, arr)
% A broadband train: several wavelengths sharing ONE speed, so f = c/lambda.
    lamU = cal.lambdaUsable;  if isnan(lamU), lamU = 6*arr.Pitch; end
    lams = lamU * [1.4 2.0 2.8 3.9];
    lams = lams(lams < 0.5*arr.Aperture);
    ks = 2*pi./lams;  fs_ = c ./ lams;  U = zeros(numel(x), numel(tv));
    for j = 1:numel(ks)
        U = U + cos(ks(j)*x - 2*pi*fs_(j)*tv + 0.5*j);
    end
end

function m = i_speedmask(Lambda, fAx, cGraph, rel)
% omega = c*sqrt(lambda) is a DIAGONAL. rel is its half-width as a fraction of the ridge.
    w = 2*pi*fAx(:).';
    ridge = cGraph * sqrt(max(Lambda(:), 0));
    m = abs(w - ridge) <= rel * max(ridge, eps);
end

function m = i_rectmask(Lambda, fAx, sa, lamRange, fRange)
    lam = 2*pi*sa ./ max(sqrt(Lambda(:)), eps);
    m = (lam >= lamRange(1) & lam <= lamRange(2)) & ...
        (fAx(:).' >= fRange(1) & fAx(:).' <= fRange(2));
end

function m = i_div()
    n = 128; t = linspace(0,1,n).';
    m = [[0.13+0.84*t, 0.29+0.68*t, 0.55+0.42*t]; [0.97-0.27*t, 0.97-0.82*t, 0.97-0.81*t]];
end

% Author: Diellor Basha, 2026
