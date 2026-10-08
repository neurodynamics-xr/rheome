function out = sensor_wavelets(outDir, varargin)
% DEMOS.SENSOR_WAVELETS  Scale, rate and speed: three controls, for GENERATION and ANALYSIS.
%
%   out = rheome.demos.sensor_wavelets()
%   out = rheome.demos.sensor_wavelets('_figures')       % also export the 5 PNGs
%
% The joint bank composes three factor lists -- W = K(lambda,omega) * psi_g(lambda) *
% psi_t(omega) -- and those three factors ARE the three controls:
%
%   psi_g   the graph member   -> SPATIAL SCALE, millimetres of sensor separation
%   psi_t   the time member    -> RATE, hertz
%   K       the joint kernel   -> SPEED, metres per second across the array
%
% Each is planted by generating a pattern and then RECOVERED by analysing it. Generation
% is a delta pushed through a member (rheome.filters.impulse); analysis is the ordinary readout.
% A control that cannot be recovered is decoration, so every one is measured here.
%
% ⚠ SPEED IS NOT MEASURED THE SAME WAY AS THE OTHER TWO. dispersion states that the
% members are not involved -- the ridge is fitted on the FULL (lambda, omega) energy --
% and that the retained band must span the diagonal. Fitting speed inside one graph member
% returns a confident wrong number: R2 stays above 0.98 while the slope is out by a factor
% of two to four, and the error is the SAME multiple at every speed, which is what marks
% it as a conditioning artefact rather than noise. So scale and rate are planted per
% member; speed is planted with a joint kernel across all wavenumbers.
%
% ⚠ NO SPEED WITHOUT rheome.sensors.calibrate. The Laplacian's lambda is dimensionless. Every
% millimetre and every metre per second below is a graph quantity times sqrt(alpha).
%
% OUTPUT:
%   out  .ok .checks .scale .rate .speed .independence .bank
%
% See also: rheome.sensors.generate, rheome.filters.impulse, rheome.jointfilterbank, rheome.sensors.calibrate
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    pp = inputParser;
    pp.addParameter('PanelArrays', '');     % JSON of {name, dim, x, y} in millimetres
    pp.parse(varargin{:});
    panelSrc = pp.Results.PanelArrays;
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});

    % ---------- the array, and the constant that carries metres ----------
    w = warning('off', 'sensors:calibrate:narrowFitWindow');
    cleanup = onCleanup(@() warning(w));
    arr = rheome.sensors.grid('Size', [24 24], 'Pitch', 10e-3);
    G   = rheome.sensors.graph(arr);
    M   = rheome.sensors.modes(G);
    cal = rheome.sensors.calibrate(G);
    sa  = sqrt(cal.alpha);                       % graph units -> metres
    fprintf('\n=== rheome.demos.sensor_wavelets : %s, alpha = %.4g m^2 ===\n', arr.Name, cal.alpha);

    fs = 500;  N = 512;
    f  = (0:N-1) * (fs/N);            fHalf = f(1:floor(N/2)+1);
    seed  = round(G.nV/2);
    basis = struct('Phi', M.Phi, 'Lambda', M.Lambda, 'Mass', speye(G.nV), 'nV', G.nV);
    kk    = sqrt(M.Lambda);

    gfb = rheome.graphfilterbank(M.Lambda, 'Wavelet','itersine', 'NumFilters',6, ...
                          'Transform', rheome.graphtransform.eigen(M.Phi, speye(G.nV), M.Lambda));
    Hs  = graphfilters(gfb);                     % [K x 6]
    kc  = centerWavenumbers(gfb);
    lamMM = 2*pi*sa ./ max(kc, eps) * 1e3;       % the scale ladder, millimetres
    bnd = framebounds(gfb);
    jfb = rheome.jointfilterbank(gfb, 'SignalLength', N, 'SamplingFrequency', fs);

    bands  = [4 8; 8 14; 14 24; 24 40];
    speeds = [0.15 0.30 0.60];                   % m/s across the array
    mScale = 2:6;                                % the members with a usable wavelength

    % ---------- Fig 1: the three controls, drawn on their own axes ----------
    f1 = fbd_fig(doExport, [50 50 1350 400]);
    subplot(1,3,1);
    plot(sqrt(gfb.Lambda_), Hs, 'LineWidth', 1); grid on; xlabel('k (graph)');
    ylabel('gain'); title(sprintf('psi_g: %d scales', gfb.NumMembers));
    subplot(1,3,2);
    hold on;
    for b = 1:size(bands,1)
        h = rheome.jointfilterbank.band(bands(b,1), bands(b,2));
        plot(fHalf, h(2*pi*fHalf), 'LineWidth', 1);
    end
    grid on; xlim([0 50]); xlabel('f (Hz)'); ylabel('gain'); title('psi_t: 4 rates');
    subplot(1,3,3);
    hold on;
    for i = 1:numel(speeds)
        plot(kk, speeds(i)/sa*kk/(2*pi), 'LineWidth', 1);
    end
    grid on; xlabel('k (graph)'); ylabel('f (Hz)'); ylim([0 60]);
    legend(compose('%.2f m/s', speeds(:)), 'Location','northwest');
    title('K: speed = a diagonal');
    sgtitle('sensor\_wavelets -- Fig 1: three factor lists, three controls');
    fbd_save(f1, outDir, 'sensor_wavelets_1_controls');

    % ---------- Fig 2 + SCALE: plant a member, recover its wavelength ----------
    c0 = M.Phi(seed,:).';
    scaleErr = zeros(1, numel(mScale));  scaleRec = zeros(1, numel(mScale));
    for i = 1:numel(mScale)
        cm   = Hs(:, mScale(i)) .* c0;                       % one scale, static
        e    = abs(cm).^2;
        kbar = sum(e .* kk) / max(sum(e), eps);              % energy-weighted wavenumber
        scaleRec(i) = 2*pi*sa/kbar * 1e3;
        scaleErr(i) = scaleRec(i)/lamMM(mScale(i)) - 1;
    end
    f2 = fbd_fig(doExport, [70 70 1350 380]);
    show = [2 4 6];
    for i = 1:3
        subplot(1,4,i);
        u = M.Phi * (Hs(:,show(i)) .* c0);
        scatter(G.Vertices(:,1)*1e3, G.Vertices(:,2)*1e3, 26, u, 'filled');
        axis equal tight off; caxis(max(abs(u))*[-1 1]); colormap(gca, i_div());
        title(sprintf('g%d: %.0f mm', show(i), lamMM(show(i))));
    end
    subplot(1,4,4);
    plot(lamMM(mScale), scaleRec, 'o-', 'LineWidth',1); hold on;
    plot(lamMM(mScale), lamMM(mScale), 'k--'); grid on; axis square;
    xlabel('planted (mm)'); ylabel('recovered (mm)'); title('scale round-trip');
    sgtitle('sensor\_wavelets -- Fig 2: the graph member sets the SPATIAL SCALE');
    fbd_save(f2, outDir, 'sensor_wavelets_2_scale');

    % ---------- Fig 3 + RATE: plant a band, recover its centroid ----------
    rateErr = zeros(1, size(bands,1));  rateRec = zeros(1, size(bands,1));
    Ub = cell(1, size(bands,1));
    for b = 1:size(bands,1)
        h  = rheome.jointfilterbank.band(bands(b,1), bands(b,2));
        Gk = Hs(:,4) * h(2*pi*f);                            % scale g4 x this rate
        Ub{b} = rheome.filters.impulse(basis, seed, Gk, 'js', N);
        S  = abs(fft(sum(M.Phi' * Ub{b}, 1), N)).^2;  S = S(1:floor(N/2)+1);
        rateRec(b) = sum(fHalf .* S) / max(sum(S), eps);
        rateErr(b) = rateRec(b)/mean(bands(b,:)) - 1;
    end
    f3 = fbd_fig(doExport, [90 90 1350 380]);
    subplot(1,3,1);
    tt = (0:N-1)/fs*1e3;
    hold on; for b = 1:size(bands,1), plot(tt, Ub{b}(seed,:)/max(abs(Ub{b}(seed,:))) + 2.2*(b-1)); end
    grid on; xlim([0 300]); xlabel('t (ms)'); yticks([]); title('same scale, four rates');
    subplot(1,3,2);
    hold on;
    for b = 1:size(bands,1)
        S = abs(fft(sum(M.Phi' * Ub{b}, 1), N)).^2;  S = S(1:floor(N/2)+1);
        plot(fHalf, S/max(S), 'LineWidth', 1);
    end
    grid on; xlim([0 50]); xlabel('f (Hz)'); ylabel('normalised'); title('temporal spectra');
    subplot(1,3,3);
    plot(mean(bands,2), rateRec, 'o-', 'LineWidth',1); hold on;
    plot(mean(bands,2), mean(bands,2), 'k--'); grid on; axis square;
    xlabel('planted centre (Hz)'); ylabel('recovered (Hz)'); title('rate round-trip');
    sgtitle('sensor\_wavelets -- Fig 3: the time member sets the RATE');
    fbd_save(f3, outDir, 'sensor_wavelets_3_rate');

    % ---------- Fig 4 + SPEED: plant a diagonal, fit it back ----------
    speedErr = zeros(1, numel(speeds));  speedRec = zeros(1, numel(speeds));
    dAll = cell(1, numel(speeds));
    for i = 1:numel(speeds)
        cG = speeds(i)/sa;
        Gk = rheome.filters.travwave(M.Lambda, f, cG, 0.15*cG/(2*pi));
        U  = rheome.filters.impulse(basis, seed, Gk, 'js', N);
        Cp = i_jspec(M.Phi, U, N);
        d  = dispersion(jfb, Cp);
        dAll{i}     = d;
        speedRec(i) = d.slope * sa;
        speedErr(i) = speedRec(i)/speeds(i) - 1;
    end
    f4 = fbd_fig(doExport, [110 110 1350 380]);
    subplot(1,3,1);
    hold on; for i = 1:numel(speeds), plot(dAll{i}.fRidge, dAll{i}.kbar, '.'); end
    grid on; xlabel('f (Hz)'); ylabel('kbar (graph)'); title('the fitted ridges');
    subplot(1,3,2);
    plot(speeds, speedRec, 'o-', 'LineWidth',1); hold on;
    plot(speeds, speeds, 'k--'); grid on; axis square;
    xlabel('planted (m/s)'); ylabel('recovered (m/s)'); title('speed round-trip');
    subplot(1,3,3);
    bar([dAll{1}.slopeR2 dAll{2}.slopeR2 dAll{3}.slopeR2]); grid on; ylim([0.9 1]);
    xticklabels(compose('%.2f', speeds(:))); xlabel('planted (m/s)'); ylabel('R^2');
    title('fit quality');
    sgtitle('sensor\_wavelets -- Fig 4: the joint kernel sets the SPEED');
    fbd_save(f4, outDir, 'sensor_wavelets_4_speed');

    % ---------- Fig 5 + ANALYSIS: a generated pattern lands in its own cell ----------
    igPlant = 4;  ibPlant = 2;
    hP = rheome.jointfilterbank.band(bands(ibPlant,1), bands(ibPlant,2));
    Up = rheome.filters.impulse(basis, seed, Hs(:,igPlant) * hP(2*pi*f), 'js', N);
    Cp = i_jspec(M.Phi, Up, N);
    E  = zeros(gfb.NumMembers, size(bands,1));
    for ig = 1:gfb.NumMembers
        for ib = 1:size(bands,1)
            hb = rheome.jointfilterbank.band(bands(ib,1), bands(ib,2));
            E(ig, ib) = sum(abs(Cp .* (Hs(:,ig) * hb(2*pi*fHalf))).^2, 'all');
        end
    end
    E = E / max(E(:));
    [~, iFlat] = max(E(:));  [igHit, ibHit] = ind2sub(size(E), iFlat);

    % Independence: plant EVERY (scale, rate) pair and check each lands in its own cell.
    nMiss = 0;  nCell = 0;
    for ig = mScale
        for ib = 1:size(bands,1)
            hb = rheome.jointfilterbank.band(bands(ib,1), bands(ib,2));
            Uc = rheome.filters.impulse(basis, seed, Hs(:,ig) * hb(2*pi*f), 'js', N);
            Cc = i_jspec(M.Phi, Uc, N);
            Ec = zeros(gfb.NumMembers, size(bands,1));
            for jg = 1:gfb.NumMembers
                for jb = 1:size(bands,1)
                    hj = rheome.jointfilterbank.band(bands(jb,1), bands(jb,2));
                    Ec(jg, jb) = sum(abs(Cc .* (Hs(:,jg) * hj(2*pi*fHalf))).^2, 'all');
                end
            end
            [~, k] = max(Ec(:));  [gHit, bHit] = ind2sub(size(Ec), k);
            nCell = nCell + 1;
            nMiss = nMiss + double(gHit ~= ig || bHit ~= ib);
        end
    end
    indep = nMiss / max(nCell, 1);
    f5 = fbd_fig(doExport, [130 130 1350 400]);
    subplot(1,3,1);
    scatter(G.Vertices(:,1)*1e3, G.Vertices(:,2)*1e3, 26, Up(:,32), 'filled');
    axis equal tight off; caxis(max(abs(Up(:,32)))*[-1 1]); colormap(gca, i_div());
    title(sprintf('planted: g%d, %d-%d Hz', igPlant, bands(ibPlant,1), bands(ibPlant,2)));
    subplot(1,3,2);
    imagesc(E'); axis xy; colorbar; colormap(gca, parula);
    xlabel('scale member'); ylabel('rate band'); yticks(1:size(bands,1));
    yticklabels(compose('%d-%d', bands(:,1), bands(:,2)));
    title(sprintf('energy: peak at g%d, band %d', igHit, ibHit));
    subplot(1,3,3);
    plot(tt, Up(seed,:)/max(abs(Up(seed,:))), 'LineWidth',1); grid on;
    xlim([0 300]); xlabel('t (ms)'); title('at the seed');
    sgtitle('sensor\_wavelets -- Fig 5: ANALYSIS puts the pattern back in its own cell');
    fbd_save(f5, outDir, 'sensor_wavelets_5_analysis');

    % ---------- self-checks ----------
    chk(end+1) = fbd_check('bank is a frame (A>0)', double(bnd.A > 0), 1, 0.5, '');
    for i = 1:numel(mScale)
        chk(end+1) = fbd_check(sprintf('scale g%d recovered', mScale(i)), ...
            scaleErr(i), 0, 0.08, 'rel'); %#ok<AGROW>
    end
    for b = 1:size(bands,1)
        chk(end+1) = fbd_check(sprintf('rate %d-%d Hz recovered', bands(b,1), bands(b,2)), ...
            rateErr(b), 0, 0.03, 'rel'); %#ok<AGROW>
    end
    for i = 1:numel(speeds)
        chk(end+1) = fbd_check(sprintf('speed %.2f m/s recovered', speeds(i)), ...
            speedErr(i), 0, 0.10, 'rel'); %#ok<AGROW>
    end
    chk(end+1) = fbd_check('analysis peaks at the planted scale', double(igHit == igPlant), 1, 0.5, '');
    chk(end+1) = fbd_check('analysis peaks at the planted rate',  double(ibHit == ibPlant), 1, 0.5, '');
    % ⚠ INDEPENDENCE IS A SWEEP, NOT ONE NUMBER. Reading the same speed back at one scale
    % proves nothing about the other two controls; a control pair is independent only if
    % EVERY combination lands in its own cell. This is the fraction that does not.
    chk(end+1) = fbd_check('every scale x rate cell is recovered', indep, 0, 0.001, 'frac');

    [ok, ~] = fbd_report(chk);
    if doExport
        E = struct('G',G, 'M',M, 'basis',basis, 'sa',sa, 'Hs',Hs, 'mScale',mScale, ...
                   'lamMM',lamMM, 'scaleRec',scaleRec, 'bands',bands, 'rateRec',rateRec, ...
                   'speeds',speeds, 'speedRec',speedRec, 'f',f, 'fs',fs, 'N',N, ...
                   'panelSrc',panelSrc);
        i_writejson(outDir, E);
    end
    out = struct('ok', ok, 'checks', chk, 'bank', lamMM, ...
                 'scale', struct('planted', lamMM(mScale), 'recovered', scaleRec, 'err', scaleErr), ...
                 'rate',  struct('planted', mean(bands,2).', 'recovered', rateRec, 'err', rateErr), ...
                 'speed', struct('planted', speeds, 'recovered', speedRec, 'err', speedErr), ...
                 'independence', indep);
end

function i_writejson(outDir, E)
% The three controls, for an interactive rendering outside MATLAB.
%
% ⚠ EXPORT THE PRODUCTS, NOT THE INGREDIENTS. Shipping the eigenbasis and letting a page
% rebuild the atoms would put the reconstruction outside the test suite.
%
% ⚠ NORMALISE PER FRAME, NOT PER CLIP. A front decays as it spreads, so dividing a clip by
% its global peak rounds late frames to zero and the panel renders a blank array -- which
% reads as a broken control rather than as a decaying wave.
%
% Three seeds, not one: a graph has NO TRANSLATION INVARIANCE, so the same member at the
% centre, at an edge and at a corner is three different atoms. That is the property the
% seed control exists to show, and it is why a seed is required rather than defaulted.
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    G = E.G;  M = E.M;  basis = E.basis;  Hs = E.Hs;  f = E.f;  N = E.N;

    seeds = i_seeds(G);
    kRate  = round(linspace(1, round(N/2), 32));      % the rate panel's frames
    kSpeed = round(linspace(1, round(N/2), 48));      % the speed panel's frames

    nS = numel(seeds);  nM = numel(E.mScale);  nB = size(E.bands,1);  nC = numel(E.speeds);

    scaleAtoms = zeros(nS, nM, G.nV);
    rateTrace  = zeros(nS, nB, numel(kRate));
    rateField  = zeros(nS, nB, G.nV, numel(kRate));
    speedField = zeros(nS, nC, G.nV, numel(kSpeed));

    for a = 1:nS
        v  = seeds(a).index;
        c0 = M.Phi(v,:).';
        for m = 1:nM
            scaleAtoms(a,m,:) = i_unit(M.Phi * (Hs(:, E.mScale(m)) .* c0));
        end
        for b = 1:nB
            h  = rheome.jointfilterbank.band(E.bands(b,1), E.bands(b,2));
            Ub = rheome.filters.impulse(basis, v, Hs(:,4) * h(2*pi*f), 'js', N);
            rateTrace(a,b,:)   = i_unit(Ub(v, kRate));
            rateField(a,b,:,:) = i_perframe(Ub(:, kRate));
        end
        for i = 1:nC
            cG = E.speeds(i)/E.sa;
            Us = rheome.filters.impulse(basis, v, ...
                     rheome.filters.travwave(M.Lambda, f, cG, 0.15*cG/(2*pi)), 'js', N);
            speedField(a,i,:,:) = i_perframe(Us(:, kSpeed));
        end
    end

    S = struct('fs', E.fs, 'nV', G.nV, ...
        'x', round(G.Vertices(:,1).'*1e3, 2), 'y', round(G.Vertices(:,2).'*1e3, 2), ...
        'tRate',  round((kRate-1)/E.fs,  4), ...
        'tSpeed', round((kSpeed-1)/E.fs, 4), ...
        'seeds', struct('name', {{seeds.name}}, 'index', [seeds.index]), ...
        'scale', struct('planted', round(E.lamMM(E.mScale),1), ...
                        'recovered', round(E.scaleRec,1), ...
                        'label', {compose('g%d', E.mScale(:)).'}, 'atoms', scaleAtoms), ...
        'rate',  struct('bands', E.bands, 'planted', round(mean(E.bands,2).',2), ...
                        'recovered', round(E.rateRec,2), ...
                        'traces', rateTrace, 'fields', rateField), ...
        'speed', struct('planted', E.speeds, 'recovered', round(E.speedRec,4), ...
                        'fields', speedField), ...
        'panel', i_panel(E));

    fid = fopen(fullfile(outDir, 'sensor_wavelets.json'), 'w');
    fwrite(fid, jsonencode(S));  fclose(fid);
    d = dir(fullfile(outDir, 'sensor_wavelets.json'));
    fprintf('  wrote %s (%.0f KB)\n', fullfile(outDir, 'sensor_wavelets.json'), d.bytes/1024);
end

function s = i_seeds(G)
% Centre, edge and corner. The three places a graph wavelet behaves differently, ordered
% by how much boundary the atom can feel.
    V = G.Vertices;  c = mean(V, 1);
    d = vecnorm(V - c, 2, 2);
    [~, iC] = min(d);                                   % centre
    [~, iK] = max(d);                                   % corner: furthest from centre
    % Edge: furthest from centre ALONG one axis while near the centre on the other, which
    % is a side rather than a corner.
    ax = V(:,1) - c(1);  ay = abs(V(:,2) - c(2));
    score = abs(ax) - 3*ay;
    [~, iE] = max(score);
    s = struct('name', {'centre','edge','corner'}, 'index', {iC, iE, iK});
end

function u = i_unit(x)
% One clip, one peak. For a static atom or a single trace this is the whole normalisation.
    m = max(abs(x(:)));
    if m <= 0, m = 1; end
    u = round(x / m, 2);
end

function Y = i_perframe(X)
% Each frame to its own peak, so a decaying front stays visible for the whole sweep.
    Y = zeros(size(X));
    for j = 1:size(X,2)
        m = max(abs(X(:,j)));
        if m > 0, Y(:,j) = X(:,j) / m; end
    end
    Y = round(Y, 2);
end

function P = i_panel(E)
% Every generator, over the SAME two controls, for the arrays the report's top panel draws.
%
% ⭐ SCALE AND RATE ARE NOT THE WAVELET'S CONTROLS. THEY ARE EVERY GENERATOR'S. sen_pattern
% has always taken a wavelength and a frequency; the report simply froze them at one value
% per array, which made a band-limited atom look like a seventh kind of pattern rather than
% the same two knobs turned differently. Exporting the grid puts them where they belong.
%
% ⚠ THE WAVELET IS AN EXCITATION, NOT A KERNEL. A travelling wave and a wavelet atom differ
% in WHAT THE KERNEL ACTS ON: the wave is an analytic field over every sensor, the atom is a
% delta at ONE sensor pushed through a band. Both carry a wavelength and a rate. That is why
% 'excitation' is a separate control here and 'wavelet' is no longer a generator.
%
% ⚠ ROTATION, SPIRAL AND TARGET ARE NOT JOINT-DOMAIN KERNELS as written. They are built from
% the polar angle about the array centre, which is a COORDINATE construct with no spectral
% counterpart -- which is exactly why they are undefined on a chain. Travelling, standing
% and colliding are functions of k.x and do have one.
%
% ⚠ 2D arrays keep all six generators; a chain keeps the three that do not need an angle.
    P = struct('arrays', {{}});
    cand = i_candidates(E.panelSrc);
    names = {'travelling','standing','rotating','spiral','target','colliding'};
    % ⚠ WHOLE CYCLES IN THE WINDOW, or the clip cannot loop. 8, 12 and 20 Hz are 2, 3 and 5
    % cycles of the 0.25 s window, so playback wraps seamlessly while still animating at
    % three visibly different speeds. An arbitrary rate leaves a jump at the wrap.
    rates = [8 12 20];                                   % Hz
    % A FIXED WINDOW, NOT A FIXED NUMBER OF CYCLES. Giving every rate its own window makes
    % each clip one cycle long, so all three animate at the same apparent speed and the
    % rate control shows nothing. One shared window is what makes a fast pattern look fast.
    Twin = 0.25;  nT = 32;                               % s, and frames: 7.8 ms apart
    tv0  = (0:nT-1) * (Twin/nT);
    for i = 1:numel(cand)
        arr = cand{i};
        if isempty(arr), continue; end
        pl  = sen_plane(arr);
        Gi  = rheome.sensors.graph(arr);  Mi = rheome.sensors.modes(Gi);
        cal = rheome.sensors.calibrate(Gi);
        lam = i_scaleladder(arr, cal);                   % metres
        bi  = struct('Phi',Mi.Phi, 'Lambda',Mi.Lambda, 'Mass',speye(Gi.nV), 'nV',Gi.nV);
        gi  = rheome.graphfilterbank(Mi.Lambda, 'Wavelet','itersine', 'NumFilters',6, ...
                  'Transform', rheome.graphtransform.eigen(Mi.Phi, speye(Gi.nV), Mi.Lambda));
        Hi  = graphfilters(gi);
        sai = sqrt(cal.alpha);
        kc  = centerWavenumbers(gi);
        lamBank = 2*pi*sai ./ max(kc, eps);              % the bank's own ladder, metres
        c = mean(Gi.Vertices, 1);
        [~, v] = min(vecnorm(Gi.Vertices - c, 2, 2));

        keep = {};  fields = {};
        for g = 1:numel(names)
            F = zeros(numel(lam), numel(rates), Gi.nV, nT);
            ok = true;
            for a = 1:numel(lam)
                for b = 1:numel(rates)
                    [U, def] = sen_pattern(arr, pl, names{g}, lam(a), rates(b), tv0);
                    if ~def, ok = false; break; end
                    F(a,b,:,:) = i_perframe(U);
                end
                if ~ok, break; end
            end
            if ok, keep{end+1} = names{g};  fields{end+1} = F; end %#ok<AGROW>
        end

        % The delta excitation, on the same grid. Its scale is the bank member nearest the
        % wavelength asked for, so the two excitations are labelled with ONE ladder.
        A = zeros(numel(lam), numel(rates), Gi.nV, nT);
        for a = 1:numel(lam)
            [~, m] = min(abs(lamBank - lam(a)));
            for b = 1:numel(rates)
                fAx = (0:nT-1) / Twin;                   % the window's own bin grid
                h = rheome.jointfilterbank.band(0.7*rates(b), 1.4*rates(b));
                U = rheome.filters.impulse(bi, v, Hi(:,m) * h(2*pi*fAx), 'js', nT);
                A(a,b,:,:) = i_perframe(U);
            end
        end

        P.arrays{end+1} = struct('name', arr.Name, 'nV', Gi.nV, 'seed', v, ...
            'scales', round(lam*1e3, 1), 'rates', rates, 'nT', nT, ...
            't', round(tv0, 4), ...
            'generators', {keep}, 'fields', {fields}, 'atom', A);
        fprintf('  panel: %-14s nV=%3d  %d generators  scales %s mm\n', ...
                arr.Name, Gi.nV, numel(keep), mat2str(round(lam*1e3,1)));
    end
end

function lam = i_scaleladder(arr, cal)
% Three wavelengths the array can actually carry: from just above the resolution floor to
% a good fraction of the aperture. Outside that a "scale control" is showing aliasing at
% one end and a single cycle at the other.
    % ⚠ ANCHOR ON PITCH AND APERTURE, NOT ON lambdaUsable ALONE. Where the usable floor sits
    % close to the aperture the two collapse onto each other and the ladder spans a factor
    % of 1.1 -- three labels that are one scale, which is not a control.
    hi = 0.75 * arr.Aperture;
    lo = max(3 * arr.Pitch, 0.25 * hi);
    lam = exp(linspace(log(hi), log(lo), 3));            % coarse -> fine
end

function cand = i_candidates(src)
% Either the coordinates the report already draws, or this machine's own constructors.
    if ~isempty(src) && exist(src, 'file')
        raw = jsondecode(fileread(src));
        cand = cell(1, numel(raw));
        for i = 1:numel(raw)
            r = raw(i);
            Pxy = [r.x(:), r.y(:)] * 1e-3;              % the report carries millimetres
            cand{i} = rheome.sensors.positions(r.name, Pxy, 'Dim', r.dim, 'Kind', r.kind);
        end
        return;
    end
    cand = {rheome.sensors.grid('utah'), rheome.sensors.grid('ecog'), i_megarray()};
end

function arr = i_megarray()
    arr = [];
    try
        arr = rheome.sensors.meg('subject01');
    catch
        fprintf('  (subject01 not cached -- MEG comes from the report coordinates or not at all)\n');
    end
end

function Cp = i_jspec(Phi, U, N)
% The joint spectrum of a field: modes along one axis, DFT bins along the other. The
% positive half only, which is the grid a 'positive' bank reports.
    Cf = fft(Phi' * U, N, 2);
    Cp = Cf(:, 1:floor(N/2)+1);
end

function cm = i_div()
% A balanced diverging map. Zero is the middle, so a sign change reads as one.
    n = 256;  x = linspace(-1, 1, n).';
    r = min(1, max(0, 1 + x));  b = min(1, max(0, 1 - x));  g = 1 - abs(x);
    cm = [r, g*0.85 + 0.15, b];
end

% Author: Diellor Basha, 2026
