function out = filterbank_joint(outDir)
% DEMOS.FILTERBANK_JOINT  jointfilterbank, in the same five figures as the other two.
%
%   rheome.demos.filterbank_joint()
%   out = rheome.demos.filterbank_joint(outDir)     % also export figures as PNGs
%
%   1 design      Bands / JointKernels / speedkernel -- what the joint parameters change
%   2 properties  centres on both axes, frame bounds, and what materialising would cost
%   3 field       a rigidly rotating pattern and a standing oscillation
%   4 decompose   joint spectrum, scalogram, energy per bin
%   5 recovery    fitted speed vs planted, and the clipped-band failure mode
%
% ⭐ GROUND TRUTH. Sectoral harmonics of degree l rotating rigidly at Omega put energy at
% temporal frequency omega_l = l*Omega and wavenumber k_l = sqrt(l(l+1))/R, so the locus
% of those points has a slope known BEFORE the fit runs. Measured here: -1.7% error with
% R^2 = 1.0000.
%
% ⚠ THE BAND MUST SPAN THE DIAGONAL. Fig 5 runs the SAME field through a spanning and a
% deliberately clipped frequency retention, so the failure mode .fSupport exists to warn
% about is visible rather than theoretical.
%
% ⭐ VOCABULARY. Both axes appear here and each takes its own correct word:
% centerWavenumbers is rad/m on the graph axis, centerFrequencies is Hz on the time axis.
% "Frequency" was wrong for lambda; on omega it is exactly right.
%
% See also: rheome.demos.filterbank_cwt, rheome.demos.filterbank_graph, rheome.jointfilterbank
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);

    % ico4/K=100 gives 4.5 vertices per half-wavelength, so rheome.eigen.modes raises no
    % resolution warning and degrees 3..7 sit well inside the trustworthy part.
    S  = fbd_sphere(4, 100);
    N  = 256;  fs = 64;  tv = (0:N-1)/fs;
    Om = 2*pi*3;                                   % 3 revolutions per second
    degs = 3:7;
    fprintf('\n=== rheome.demos.filterbank_joint : %d verts, K = %d, %d samples at %g Hz ===\n', ...
        S.nV, numel(S.Lambda), N, fs);

    gfb = rheome.graphfilterbank(S.Lambda, 'Wavelet','itersine', 'NumFilters',6, 'Transform', S.T);

    % ---------- Fig 1: what the joint parameters do ----------
    f1 = fbd_fig(doExport, [60 60 1250 400]);
    subplot(1,3,1);
    jb = rheome.jointfilterbank(gfb, 'SignalLength',N, 'SamplingFrequency',fs, 'Bands',[2 8; 10 20]);
    plot(centerFrequencies(jb), '.-', 'MarkerSize', 12); grid on;
    xlabel('member #'); ylabel('centre frequency (Hz)');
    title(sprintf('Bands: %d graph \\times %d time = %d members', ...
        jb.NumGraph, jb.NumTime, jb.NumMembers));
    subplot(1,3,2);
    js = rheome.jointfilterbank(gfb, 'SignalLength',N, 'SamplingFrequency',fs, ...
                         'JointKernels', {rheome.jointfilterbank.speedkernel(2, 15)});
    spectralResponse(js, 3);
    title('speedkernel: a DIAGONAL member');
    subplot(1,3,3);
    jn = rheome.jointfilterbank(gfb, 'SignalLength',N, 'SamplingFrequency',fs);
    spectralResponse(jn, 3);
    title('separable member: a HORIZONTAL band');
    sgtitle('rheome.jointfilterbank -- Fig 1: design-stage parameters  (cf. cwt Fig 1)');
    fbd_save(f1, outDir, 'joint_1_design');

    % ---------- Fig 2: properties ----------
    jfb = rheome.jointfilterbank(gfb, 'SignalLength',N, 'SamplingFrequency',fs);
    bj  = framebounds(jfb);
    bytesIfMaterialised = numel(jfb.Lambda) * jfb.NumOmega * jfb.NumMembers * 8;

    f2 = fbd_fig(doExport, [70 70 1250 420]);
    subplot(1,3,1); plot(centerWavenumbers(jfb), '.-', 'MarkerSize', 12); grid on;
    xlabel('member #'); ylabel('k (rad/m)'); title('centerWavenumbers -- rad/m');
    subplot(1,3,2); plot(centerFrequencies(jfb), '.-', 'MarkerSize', 12); grid on;
    xlabel('member #'); ylabel('f (Hz)'); title('centerFrequencies -- Hz IS correct here');
    subplot(1,3,3);
    i_kimage(jfb.Frequencies, sqrt(jfb.Lambda), bj.S); colorbar;
    title(sprintf('frame operator  A=%.3g  B=%.3g', bj.A, bj.B));
    sgtitle(sprintf(['rheome.jointfilterbank -- Fig 2: properties. Materialised the bank would be ' ...
        '%.0f kB; it is stored as FACTORS instead'], bytesIfMaterialised/1e3));
    fbd_save(f2, outDir, 'joint_2_properties');

    % ---------- Fig 3: spatiotemporal fields with known features ----------
    [U, gt] = fbd_rotating(S, degs, Om, tv);
    lStand = 5;  fStand = 12;
    Ustand = fbd_sectoral(S, lStand, 0) * cos(2*pi*fStand*tv);
    surfS = struct('Vertices', S.V, 'Faces', S.F, 'nV', S.nV);

    f3 = fbd_fig(doExport, [80 80 1250 560]);
    snaps = [1 6 11];
    for i = 1:3
        a = subplot(2,3,i);
        rheome.show.surface(surfS, U(:,snaps(i)), 'Parent', a, 'Colorbar', false);
        title(a, sprintf('rotating: t = %.3f s', tv(snaps(i))));
    end
    subplot(2,3,[4 5 6]);
    plot(tv, U(1,:), '-', tv, Ustand(1,:), '--', 'LineWidth', 1.1); grid on;
    xlabel('time (s)'); ylabel('field at vertex 1');
    legend({'rotating','standing'}, 'Location','best');
    sgtitle(sprintf(['rheome.jointfilterbank -- Fig 3: rigid rotation at %.1f rev/s, degrees %s ' ...
        '(cf. cwt Fig 3)'], Om/(2*pi), mat2str(degs)));
    fbd_save(f3, outDir, 'joint_3_field');

    % ---------- Fig 4: decomposition ----------
    C  = i_joint(S, U, N);
    Cs = i_joint(S, Ustand, N);
    Ejoint = scalogram(jfb, C);

    f4 = fbd_fig(doExport, [90 90 1250 720]);
    kAxis = sqrt(jfb.Lambda);
    subplot(2,2,1);
    i_kimage(jfb.Frequencies, kAxis, abs(C)); hold on; colorbar;
    % ⚠ The overlay MUST be mapped to rows the same way the image is, or it lands
    % somewhere else: sqrt(lambda) is non-uniform and the image is drawn by index.
    plot(gt.omega/(2*pi), i_krow(kAxis, gt.k), 'r+', 'MarkerSize', 11, 'LineWidth', 1.5);
    title('rotating: a DIAGONAL ridge  (+ = planted (k_l, f_l))');
    subplot(2,2,2);
    i_kimage(jfb.Frequencies, kAxis, abs(Cs)); colorbar;
    title('standing: a SPOT, not a diagonal');
    subplot(2,2,3);
    imagesc(jfb.Frequencies, 1:jfb.NumMembers, Ejoint); set(gca,'YDir','normal'); colorbar;
    xlabel('f (Hz)'); ylabel('member #'); title('scalogram [N_f \times n_\Omega]');
    subplot(2,2,4);
    plot(jfb.Frequencies, sum(Ejoint,1), '.-'); grid on;
    xlabel('f (Hz)'); ylabel('energy'); title('energy per frequency bin');
    sgtitle('rheome.jointfilterbank -- Fig 4: decomposition  (cf. cwt Fig 4)');
    fbd_save(f4, outDir, 'joint_4_decomposition');

    % ---------- Fig 5: recovery, and the clipped-band failure mode ----------
    d  = dispersion(jfb, C);
    % Derive the retained bins ONCE, so the clipped bank and the sliced spectrum cannot
    % drift apart -- exactly the mismatch iscompatible exists to catch.
    keep  = jfb.Frequencies >= 9 & jfb.Frequencies <= 14;
    jClip = rheome.jointfilterbank(gfb, 'Frequencies', jfb.Frequencies(keep));
    dClip = dispersion(jClip, C(:, keep));

    checks = [ ...
        fbd_check('rotating: fitted speed', d.slope, gt.slope, 0.10*gt.slope, 'm/s'), ...
        fbd_check('rotating: fit quality R2', d.slopeR2, 1, 0.10, ''), ...
        fbd_check('standing: peak frequency', i_peakfreq(Cs, jfb), fStand, 1.5, 'Hz'), ...
        fbd_check('clipped support is narrower', ...
            double(diff(dClip.fSupport) < diff(d.fSupport)), 1, 0.5, '')];

    f5 = fbd_fig(doExport, [100 100 1150 440]);
    subplot(1,2,1);
    plot(d.kbar, 2*pi*d.fRidge, '.', 'MarkerSize', 9); hold on;
    plot(gt.k, gt.omega, 'r+', 'MarkerSize', 13, 'LineWidth', 1.7);
    kk = linspace(min(d.kbar), max(d.kbar), 10);
    plot(kk, gt.slope*kk, 'k--', 'LineWidth', 1.2);
    grid on; xlabel('mean wavenumber (rad/m)'); ylabel('\omega (rad/s)');
    legend({'measured ridge','planted (k_l,\omega_l)','planted slope'}, 'Location','northwest');
    title(sprintf('fitted %.4f vs planted %.4f m/s  (R^2 %.4f)', ...
        d.slope, gt.slope, d.slopeR2));
    subplot(1,2,2);
    bar([diff(d.fSupport), diff(dClip.fSupport)]); grid on;
    set(gca, 'XTickLabel', {'full axis','clipped 9-14 Hz'});
    ylabel('.fSupport width (Hz)');
    title(sprintf('clipped fit gives %.4f m/s -- a band that does not\nspan the diagonal cannot measure speed', ...
        dClip.slope));
    sgtitle('rheome.jointfilterbank -- Fig 5: measured vs analytic  (cf. cwt Fig 5)');
    fbd_save(f5, outDir, 'joint_5_recovery');

    [ok, ~] = fbd_report(checks);
    out = struct('checks', checks, 'ok', ok, ...
                 'plantedSlope', gt.slope, 'recoveredSlope', d.slope, ...
                 'slopeR2', d.slopeR2, 'fSupport', d.fSupport, ...
                 'clippedSlope', dClip.slope, 'clippedSupport', dClip.fSupport, ...
                 'standingFreqHz', i_peakfreq(Cs, jfb));
end

function C = i_joint(S, U, N)
% Joint spectrum: project onto the eigenbasis, then FFT along time, keeping the positive
% half so the axes match a 'positive' jointfilterbank.
    c    = S.T.forward(U);
    Chat = fft(c, N, 2) / sqrt(N);
    C    = Chat(:, 1:floor(N/2)+1);
end

function i_kimage(f, k, M)
% Draw a (k, f) image with rows by INDEX and wavenumbers as tick LABELS.
%
% ⚠ imagesc assumes a UNIFORMLY spaced y vector. sqrt(lambda) is not: on an ico4/K=100
% sphere the spacing varies by orders of magnitude, and passing it directly put the
% degree-3 harmonic at k = 8.8 rad/m instead of 34.6 -- an error of 25.9 rad/m that looked
% entirely plausible.
    imagesc(f, 1:numel(k), M);
    set(gca, 'YDir', 'normal');
    yt = round(linspace(1, numel(k), 6));
    set(gca, 'YTick', yt, 'YTickLabel', compose('%.0f', k(yt)));
    xlabel('f (Hz)'); ylabel('k = \surd\lambda (rad/m)');
end

function r = i_krow(k, kWanted)
% Map a wavenumber to its row on the image above.
%
% ⚠ NOT interp1: a sphere has DEGENERATE MULTIPLETS (degree l has multiplicity 2l+1), so
% sqrt(lambda) repeats and interp1 refuses non-unique sample points. Take the centre of
% the matching multiplet, which is where the energy visually sits anyway.
    k = double(k(:));
    r = zeros(size(kWanted));
    for i = 1:numel(kWanted)
        d = abs(k - double(kWanted(i)));
        r(i) = mean(find(d <= min(d) + 1e-9*max(1, max(k))));
    end
end

function f = i_peakfreq(C, jfb)
    [~, i] = max(sum(abs(C).^2, 1));
    f = jfb.Frequencies(i);
end

% Author: Diellor Basha, 2026
