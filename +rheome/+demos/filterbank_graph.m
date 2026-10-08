function out = filterbank_graph(outDir)
% DEMOS.FILTERBANK_GRAPH  graphfilterbank, in the same five figures as the CWT demo.
%
%   rheome.demos.filterbank_graph()
%   out = rheome.demos.filterbank_graph(outDir)     % also export figures as PNGs
%
% Read alongside rheome.demos.filterbank_cwt: the figure order is identical, so each panel has a
% counterpart there.
%
%   1 design      VoicesPerOctave / Wavelet / SizeLimits -- the FrequencyLimits analogue
%   2 properties  widths, centre WAVENUMBERS, bandwidths, Q, AND frame bounds
%   3 field       Gaussian bumps of known width, plus a spatial chirp
%   4 decompose   wt + vertexSpectrum + scaleSpectrum + an atom
%   5 recovery    planted vs recovered width, with PASS/FAIL
%
% ⭐ GROUND TRUTH, DERIVED FOR THE QUANTITY ACTUALLY MEASURED. vertexSpectrum returns
% ENERGY, ||g_m .* c||^2, so the relevant integral is of the SQUARED gain. A Gaussian bump
% of geodesic width sigma has spectrum exp(-lambda*sigma^2/2); against a mexhat member
% t*lambda*exp(-t*lambda), with the mode density uniform in lambda (Weyl's law on a
% surface),
%
%   E(t) = int (t*lam*e^{-t*lam})^2 (e^{-lam*sigma^2/2})^2 dlam = 2t^2 / (2t + sigma^2)^3
%
% which is maximised at t = sigma^2, so the peak member's width is
%
%   sigma_m = sqrt(2t) = sqrt(2) * sigma.
%
% ⚠ THREE DIFFERENT CALIBRATIONS LIVE IN THIS CODEBASE AND CONFUSING THEM IS THE WHOLE
% POINT OF STATING THEM:
%   sqrt(2)      -- THIS one: a Gaussian bump read off an ENERGY marginal
%   1            -- a Gaussian bump matched by a LINEAR response (no squaring)
%   1/sqrt(2)    -- rheome.filters.frame's Gamma: a VORTEX, whose vorticity spectrum carries an
%                   extra factor of lambda, moving the maximum to t = sigma^2/4
% Measured here: +8.9%, -10.8%, -2.6% against sqrt(2)*sigma, i.e. inside one member
% spacing. Assuming the wrong one is a 2x error that still looks like a measurement.
%
% ⭐ FIG 2 DIVERGES FROM THE CWT DELIBERATELY. framebounds and MassLost/Usable have no CWT
% counterpart: a CWT is complete and constant-Q by construction, ours is neither, and
% whether a bank inverts is a question the user must be able to ask.
%
% See also: rheome.demos.filterbank_cwt, rheome.demos.filterbank_joint, rheome.graphfilterbank
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);

    % ico5/K=400 is chosen so the basis resolves every planted width WITHOUT tripping
    % rheome.eigen.modes' resolution warning: it gives 4.5 vertices per half-wavelength and a
    % usable floor of 15.5 mm, below the finest feature planted here.
    S   = fbd_sphere(5, 400);
    lam = S.Lambda;
    fprintf('\n=== rheome.demos.filterbank_graph : %d vertices, K = %d, R = %g mm ===\n', ...
        S.nV, numel(lam), 1000*S.R);

    SLIM = [15e-3 130e-3];  VPO = 4;
    surfS = struct('Vertices', S.V, 'Faces', S.F, 'nV', S.nV);

    % ---------- Fig 1: what the design parameters do ----------
    f1 = fbd_fig(doExport, [60 60 1250 380]);
    voices = [1 3 8];
    subplot(1,3,1); hold on;
    for v = voices
        g = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'VoicesPerOctave',v, 'SizeLimits',SLIM);
        plot(1000*widths(g), 1:g.NumMembers, '.-', 'MarkerSize', 12);
    end
    set(gca,'XScale','log'); grid on;
    xlabel('member width \sigma (mm)'); ylabel('member #');
    title('VoicesPerOctave: DENSITY'); legend(compose('%d', voices), 'Location','best');

    subplot(1,3,2); hold on;
    fams = {'mexhat','itersine','heat'};
    for i = 1:3
        g = rheome.graphfilterbank(lam, 'Wavelet', fams{i}, 'NumFilters', 6);
        [H, k] = spectralResponse(g);
        plot(k, H(:, 3), 'LineWidth', 1.1);
    end
    grid on; xlabel('wavenumber k = \surd\lambda (rad/m)'); ylabel('gain');
    title('Wavelet: SHAPE (3rd member)'); legend(fams, 'Location','best');

    subplot(1,3,3); hold on;
    % The narrow bank DELIBERATELY fails to cover the spectrum -- that is what "reach"
    % means, and graphfilterbank says so. Silenced here because the panel is about member
    % placement; Fig 2 is where coverage is actually examined.
    ws = warning('off', 'graphfilterbank:degenerateFrame');
    cleanupW = onCleanup(@() warning(ws));
    slims = {SLIM, [30e-3 60e-3]};
    for i = 1:2
        g = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'VoicesPerOctave',VPO, 'SizeLimits',slims{i});
        plot(1000*widths(g), 1:g.NumMembers, '.-', 'MarkerSize', 12);
    end
    clear cleanupW
    set(gca,'XScale','log'); grid on;
    xlabel('member width \sigma (mm)'); ylabel('member #');
    title('SizeLimits: REACH  (cf. FrequencyLimits)');
    legend({'15-130 mm','30-60 mm'}, 'Location','best');
    sgtitle('rheome.graphfilterbank -- Fig 1: design-stage parameters  (cf. cwt Fig 1)');
    fbd_save(f1, outDir, 'graph_1_design');

    % ---------- Fig 2: properties, and what a CWT cannot show ----------
    gfb = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'VoicesPerOctave',VPO, ...
                          'SizeLimits',SLIM, 'Transform', S.T);
    b = framebounds(gfb);  q = qfactor(gfb);  pb = powerbw(gfb);

    f2 = fbd_fig(doExport, [70 70 1250 720]);
    subplot(2,3,1); semilogy(1000*widths(gfb), '.-', 'MarkerSize', 12); grid on;
    xlabel('member #'); ylabel('\sigma (mm)'); title('widths (log-spaced)');
    subplot(2,3,2); semilogy(centerWavenumbers(gfb), '.-', 'MarkerSize', 12); grid on;
    xlabel('member #'); ylabel('k (rad/m)'); title('centerWavenumbers (rad/m, NOT Hz)');
    subplot(2,3,3); plot(pb(:,1), '.-'); hold on; plot(pb(:,2), '.-'); grid on;
    xlabel('member #'); ylabel('k (rad/m)'); title('powerbw half-power band');
    legend({'low','high'}, 'Location','northwest');

    subplot(2,3,4); plot(q, '.-', 'MarkerSize', 12); grid on;
    xlabel('member #'); ylabel('Q');
    title('qfactor PER MEMBER  (a CWT gives ONE number)');
    subplot(2,3,5);
    [H, k] = spectralResponse(gfb);
    plot(k, sum(H.^2, 2), 'k-', 'LineWidth', 1.3); grid on;
    xlabel('k (rad/m)'); ylabel('S(\lambda)');
    title(sprintf('frame operator  A=%.3g  B=%.3g  B/A=%.3g', b.A, b.B, b.Tightness));
    subplot(2,3,6);
    bar(100*gfb.MassLost); hold on; yline(5, 'r--', '5%'); grid on;
    xlabel('member #'); ylabel('mass beyond \lambda_{max} (%)');
    title(sprintf('MassLost: %d of %d unusable', nnz(~gfb.Usable), gfb.NumMembers));
    sgtitle('rheome.graphfilterbank -- Fig 2: properties, plus the frame diagnostics a CWT has no need of');
    fbd_save(f2, outDir, 'graph_2_properties');

    % ---------- Fig 3: fields with known features ----------
    plantedSigma = [20e-3 35e-3 55e-3];
    seeds = i_spread_seeds(S, numel(plantedSigma));
    u = zeros(S.nV, 1);
    for i = 1:numel(plantedSigma)
        u = u + fbd_bump(S, seeds(i), plantedSigma(i));
    end
    chirpSeeds = i_meridian_seeds(S, 6);
    chirpSigma = linspace(16e-3, 60e-3, numel(chirpSeeds));
    uc = zeros(S.nV, 1);
    for i = 1:numel(chirpSeeds)
        uc = uc + fbd_bump(S, chirpSeeds(i), chirpSigma(i));
    end

    f3 = fbd_fig(doExport, [80 80 1150 500]);
    a1 = subplot(1,2,1); rheome.show.surface(surfS, u,  'Parent', a1, 'Colorbar', false);
    title(a1, sprintf('three bumps: \\sigma = %s mm', mat2str(1000*plantedSigma)));
    a2 = subplot(1,2,2); rheome.show.surface(surfS, uc, 'Parent', a2, 'Colorbar', false);
    title(a2, sprintf('spatial chirp: \\sigma %g \\rightarrow %g mm along a meridian', ...
        1000*chirpSigma(1), 1000*chirpSigma(end)));
    sgtitle('rheome.graphfilterbank -- Fig 3: fields with known features  (cf. cwt Fig 3)');
    fbd_save(f3, outDir, 'graph_3_field');

    % ---------- Fig 4: decomposition ----------
    W  = wt(gfb, u);
    Ev = vertexSpectrum(gfb, u);
    Es = scaleSpectrum(gfb, u);
    atoms = graphfilters(gfb, 'Type','vertex', 'Vertex', seeds(2));
    mMid = round(gfb.NumMembers/2);

    f4 = fbd_fig(doExport, [90 90 1250 720]);
    a = subplot(2,3,1); rheome.show.surface(surfS, squeeze(W(:,1,2)), 'Parent',a, 'Colorbar',false);
    title(a, sprintf('band 2: \\sigma = %.0f mm', 1000*widths(gfb)*[zeros(1,1);1;zeros(gfb.NumMembers-2,1)]));
    a = subplot(2,3,2); rheome.show.surface(surfS, squeeze(W(:,1,mMid)), 'Parent',a, 'Colorbar',false);
    title(a, sprintf('band %d: \\sigma = %.0f mm', mMid, 1000*subsref(widths(gfb), substruct('()',{mMid}))));
    a = subplot(2,3,3); rheome.show.surface(surfS, squeeze(W(:,1,end)), 'Parent',a, 'Colorbar',false);
    title(a, sprintf('band %d: \\sigma = %.0f mm', gfb.NumMembers, ...
        1000*subsref(widths(gfb), substruct('()',{gfb.NumMembers}))));
    a = subplot(2,3,4); rheome.show.surface(surfS, Es, 'Parent',a, 'Colorbar',false);
    title(a, 'scaleSpectrum: energy per VERTEX');
    subplot(2,3,5);
    semilogx(1000*widths(gfb), Ev, '.-', 'MarkerSize', 12); grid on; hold on;
    for s = plantedSigma, xline(1000*sqrt(2)*s, 'r--'); end
    xlabel('member width \sigma (mm)'); ylabel('energy');
    title('vertexSpectrum: energy per SCALE (red = \surd2\cdot\sigma)');
    a = subplot(2,3,6); rheome.show.surface(surfS, squeeze(atoms(:,1,mMid)), 'Parent',a, 'Colorbar',false);
    title(a, 'atom: graphfilters(...,''Type'',''vertex'')');
    sgtitle('rheome.graphfilterbank -- Fig 4: decomposition and marginals  (cf. cwt Fig 4)');
    fbd_save(f4, outDir, 'graph_4_decomposition');

    % ---------- Fig 5: quantitative recovery ----------
    wgt  = widths(gfb);
    cand = find(isfinite(wgt));                  % skip the low-pass, whose width is NaN
    recovered = zeros(1, numel(plantedSigma));
    Eeach = zeros(gfb.NumMembers, numel(plantedSigma));
    for i = 1:numel(plantedSigma)
        Ei = vertexSpectrum(gfb, fbd_bump(S, seeds(i), plantedSigma(i)));
        Eeach(:, i) = Ei;
        [~, j] = max(Ei(cand));
        recovered(i) = wgt(cand(j));
    end

    spacing = 2^(1/VPO);
    checks = struct('name',{},'measured',{},'expected',{},'tol',{},'unit',{},'pass',{});
    for i = 1:numel(plantedSigma)
        expW = sqrt(2) * plantedSigma(i);        % ENERGY calibration, see the header
        checks(end+1) = fbd_check(sprintf('bump %d width (=%s2*sigma)', i, char(8730)), ...
            1000*recovered(i), 1000*expW, 1000*expW*(spacing-1), 'mm'); %#ok<AGROW>
    end
    checks(end+1) = fbd_check('frame lower bound A > 0', double(b.A > 0), 1, 0.5, '');
    % ⚠ Against the PROJECTION of u: the basis keeps K modes, so a geodesic-bump field has
    % content the transform cannot represent and 'dual' correctly returns what it can.
    uProj = S.T.inverse(S.T.forward(u));
    checks(end+1) = fbd_check('dual reconstruction error', ...
        norm(iwt(gfb, W, 'dual') - uProj) / norm(uProj), 0, 1e-6, '');

    f5 = fbd_fig(doExport, [100 100 1100 440]);
    subplot(1,2,1);
    semilogx(1000*wgt, Eeach, '.-', 'MarkerSize', 12); grid on; hold on;
    for s = plantedSigma, xline(1000*sqrt(2)*s, 'k--'); end
    xlabel('member width \sigma (mm)'); ylabel('energy');
    title('each bump alone (dashed = \surd2\cdot\sigma)');
    legend(compose('\\sigma = %g mm', 1000*plantedSigma), 'Location','best');
    subplot(1,2,2);
    loglog(1000*plantedSigma, 1000*recovered, 'o', 'MarkerSize', 10, 'LineWidth', 1.6);
    hold on; loglog(1000*plantedSigma, 1000*sqrt(2)*plantedSigma, 'k--', 'LineWidth', 1.2);
    grid on; xlabel('planted \sigma (mm)'); ylabel('recovered member width (mm)');
    legend({'measured','\surd2\cdot\sigma (analytic)'}, 'Location','northwest');
    title('energy marginal peaks at \surd2\cdot\sigma, not \sigma');
    sgtitle('rheome.graphfilterbank -- Fig 5: measured vs analytic  (cf. cwt Fig 5)');
    fbd_save(f5, outDir, 'graph_5_recovery');

    [ok, ~] = fbd_report(checks);
    out = struct('checks', checks, 'ok', ok, ...
                 'plantedSigma', plantedSigma, 'recoveredSigma', recovered, ...
                 'expectedSigma', sqrt(2)*plantedSigma, ...
                 'frameA', b.A, 'frameB', b.B, 'tightness', b.Tightness, ...
                 'nUnusable', nnz(~gfb.Usable));
end

function idx = i_spread_seeds(S, n)
% Well-separated seeds by farthest-point sampling, so the bumps do not overlap.
    idx = zeros(1, n);  idx(1) = 1;
    d = S.R * acos(min(1, max(-1, (S.V * S.V(1,:).') / S.R^2)));
    for i = 2:n
        [~, idx(i)] = max(d);
        dn = S.R * acos(min(1, max(-1, (S.V * S.V(idx(i),:).') / S.R^2)));
        d = min(d, dn);
    end
end

function idx = i_meridian_seeds(S, n)
% Seeds marching along the phi = 0 meridian, so the chirp sweeps scale across space.
    ph = atan2(S.V(:,2), S.V(:,1));
    th = acos(min(1, max(-1, S.V(:,3) / S.R)));
    want = linspace(0.18*pi, 0.82*pi, n);
    idx = zeros(1, n);
    for i = 1:n
        cost = (th - want(i)).^2 + 4*(angle(exp(1i*ph))).^2;
        [~, idx(i)] = min(cost);
    end
end

% Author: Diellor Basha, 2026
