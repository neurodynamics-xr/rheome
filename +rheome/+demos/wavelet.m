function out = wavelet(SurfaceFile, outDir, K)
% DEMOS.WAVELET  Localize spectral filters as atoms on the cortical surface.
%
%   rheome.demos.wavelet(SurfaceFile)
%   rheome.demos.wavelet(SurfaceFile, outDir)     % also export figures as PNGs to outDir
%   out = rheome.demos.wavelet(SurfaceFile, outDir, K)
%
% The cortical-surface analogue of GSPBox's gsp_demo_wavelet: instead of a graph we
% use the Laplace-Beltrami eigenbasis of a real cortex, and instead of a Kronecker on
% a graph vertex we drop a delta on a cortical vertex and filter it. Two experiments,
% mirroring the reference demo:
%
%   1. HEAT DIFFUSION -- a heat/low-pass kernel localized at one vertex, at four
%      diffusion times. Each atom is a smooth positive bump that widens with t.
%   2. MEXICAN-HAT WAVELETS -- a band-pass wavelet filterbank localized at the same
%      vertex, at four scales. Each atom is a signed core-plus-surround at its own
%      spatial scale (a cortical wavelet atom).
%
% Because a whole cortex is two disconnected hemispheres, an atom seeded on one
% hemisphere stays entirely on that hemisphere -- the diffusion cannot cross the gap,
% which is exactly the anatomy-respecting behaviour we want.
%
% INPUTS:
%   SurfaceFile : path to a Brainstorm cortex .mat
%   outDir      : optional folder; if given, each figure is exported as a PNG
%   K           : number of eigenmodes (default 500). More modes -> sharper atoms.
%
% OUTPUT (optional struct out): basis, seed vertex, the two gain banks, the atoms.
%
% Uses: rheome.io.read.surface, rheome.operators.laplace_beltrami, rheome.eigen.modes,
%       rheome.filters.heat, rheome.filters.mexhat, rheome.filters.localize, rheome.show.surface
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(SurfaceFile)
        error('demos:wavelet:args', 'Pass the path to a Brainstorm cortex surface .mat.');
    end
    if nargin < 2, outDir = ''; end
    if nargin < 3 || isempty(K), K = 500; end
    doExport = ~isempty(outDir);
    if doExport && ~exist(outDir, 'dir'), mkdir(outDir); end

    fprintf('\n=== Atom-design demo: localizing filters as cortical atoms ===\n');

    % --- 1. surface, operators, eigenbasis ---
    S = rheome.io.read.surface(SurfaceFile);
    fprintf('Surface: %s  (%d vertices)\n', S.Comment, S.nV);
    [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces, 'galerkin');
    fprintf('Computing %d Laplace-Beltrami eigenmodes...\n', K);
    basis = rheome.eigen.modes(L, M, K);
    lmax  = max(basis.Lambda);
    fprintf('  spectrum: lambda in [%.3g, %.3g]\n', basis.Lambda(1), lmax);

    % --- 2. seed vertex: most-lateral vertex of the right hemisphere (Y<0 in SCS) ---
    % SCS axes: +X anterior, +Y left ear, +Z up. The right hemisphere is at Y<0; its
    % most-lateral vertex sits centre-stage in a right-lateral camera view.
    [~, seed] = min(S.Vertices(:, 2));
    camView = [0 8];   % right-lateral
    fprintf('Seed vertex: %d   (right-lateral)\n', seed);

    % --- 3. HEAT diffusion at four times (spectrum-normalized, dimensionless) ---
    heatT   = [1 4 12 30];
    gHeat   = rheome.filters.heat(basis.Lambda, heatT, lmax);      % [K x 4]
    aHeat   = rheome.filters.localize(basis, seed, gHeat);         % [nV x 4]

    fH = figure('Color', 'w', 'Visible', ternary(doExport,'off','on'), 'Position', [80 80 900 780]);
    for ii = 1:4
        ax = subplot(2, 2, ii, 'Parent', fH);
        rheome.show.surface(S, aHeat(:, ii), 'Parent', ax, 'View', camView, ...
                     'Title', sprintf('Heat diffusion  t = %g', heatT(ii)));
    end
    if doExport, exportgraphics(fH, fullfile(outDir, 'wavelet_heat.png'), 'Resolution', 130); end

    % --- 4. MEXICAN-HAT wavelet filterbank at four scales ---
    % Place the band-pass peaks (lambda = 1/t) across the spectrum: low frequency
    % (coarse, wide atom) -> high frequency (fine atom).
    lamCenters = lmax * [0.05 0.15 0.40 0.80];
    mexT       = 1 ./ lamCenters;
    gMex       = rheome.filters.mexhat(basis.Lambda, mexT);        % [K x 4]
    aMex       = rheome.filters.localize(basis, seed, gMex);       % [nV x 4]

    fW = figure('Color', 'w', 'Visible', ternary(doExport,'off','on'), 'Position', [100 100 900 780]);
    for ii = 1:4
        ax = subplot(2, 2, ii, 'Parent', fW);
        rheome.show.surface(S, aMex(:, ii), 'Parent', ax, 'View', camView, ...
                     'Title', sprintf('Wavelet  scale %d  (\\lambda_c = %.0f)', ii, lamCenters(ii)));
    end
    if doExport, exportgraphics(fW, fullfile(outDir, 'wavelet_atoms.png'), 'Resolution', 130); end

    % --- 5. filterbank spectra (like gsp_plot_filter) ---
    fS = figure('Color', 'w', 'Visible', ternary(doExport,'off','on'), 'Position', [120 120 900 380]);
    ax1 = subplot(1, 2, 1, 'Parent', fS);
    plot(ax1, basis.Lambda, gHeat, 'LineWidth', 1.5);
    xlabel(ax1, '\lambda'); ylabel(ax1, 'g(\lambda)'); title(ax1, 'Heat kernels'); grid(ax1, 'on');
    legend(ax1, arrayfun(@(t) sprintf('t=%g', t), heatT, 'uni', 0), 'Location', 'northeast');
    ax2 = subplot(1, 2, 2, 'Parent', fS);
    plot(ax2, basis.Lambda, gMex, 'LineWidth', 1.5);
    xlabel(ax2, '\lambda'); ylabel(ax2, 'g(\lambda)'); title(ax2, 'Mexican-hat filterbank'); grid(ax2, 'on');
    legend(ax2, arrayfun(@(k) sprintf('scale %d', k), 1:4, 'uni', 0), 'Location', 'northeast');
    if doExport, exportgraphics(fS, fullfile(outDir, 'wavelet_spectra.png'), 'Resolution', 130); end

    fprintf('Done.%s\n\n', ternary(doExport, sprintf(' Figures exported to %s', outDir), ''));

    if nargout > 0
        out = struct('basis', basis, 'seed', seed, 'heatT', heatT, 'gHeat', gHeat, ...
                     'atomsHeat', aHeat, 'mexScales', lamCenters, 'gMex', gMex, 'atomsMex', aMex);
    end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

% Author: Diellor Basha, 2026
