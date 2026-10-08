function out = filterbank_cwt(outDir)
% DEMOS.FILTERBANK_CWT  What a filterbank's parameters do, on MATLAB's validated class.
%
%   rheome.demos.filterbank_cwt()
%   out = rheome.demos.filterbank_cwt(outDir)     % also export figures as PNGs
%
% The REFERENCE of the three-demo suite. It establishes the five-figure order that
% rheome.demos.filterbank_graph and rheome.demos.filterbank_joint then repeat on our own classes, so
% the analogy can be read by scrolling:
%
%   1 design      what VoicesPerOctave / Wavelet / FrequencyLimits actually change
%   2 properties  scales, centre frequencies, bandwidths, Q
%   3 signal      a chirp and a burst, with their known features annotated
%   4 decompose   scalogram + both marginals
%   5 recovery    measured vs analytic, with PASS/FAIL
%
% GROUND TRUTH. The chirp sweeps 20 -> 200 Hz linearly over 1 s, so its instantaneous
% frequency is known at every instant and the scalogram ridge must track it. The burst is
% a Gaussian-windowed 100 Hz tone at t = 0.5 s, so the scalogram maximum must land there.
%
% ⭐ NOTE FOR THE COMPARISON. qfactor here returns ONE number: a CWT is CONSTANT-Q by
% construction, which is the defining property of a wavelet transform. graphfilterbank
% returns Q per member, because a mexhat bank over a TRUNCATED spectrum is not constant-Q
% near the edges. Fig 2 of each demo shows the contrast.
%
% Requires the Wavelet Toolbox; without it the demo SKIPS rather than fails.
%
% See also: rheome.demos.filterbank_graph, rheome.demos.filterbank_joint, cwtfilterbank
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    out = struct('checks', [], 'ok', false, 'skipped', false);

    if isempty(ver('wavelet'))
        fprintf('rheome.demos.filterbank_cwt: Wavelet Toolbox not installed -- SKIPPED.\n');
        out.skipped = true;  out.ok = true;
        return;
    end

    fs = 1000;  N = 1024;  t = (0:N-1)/fs;
    fprintf('\n=== rheome.demos.filterbank_cwt : %d samples at %g Hz ===\n', N, fs);

    % ---------- Fig 1: what the design parameters do ----------
    f1 = fbd_fig(doExport, [60 60 1250 380]);
    voices = [4 10 20];
    subplot(1,3,1); hold on;
    for v = voices
        fb = cwtfilterbank('SignalLength',N,'SamplingFrequency',fs,'VoicesPerOctave',v);
        plot(centerFrequencies(fb), 1:numel(scales(fb)), '.-', 'MarkerSize', 8);
    end
    set(gca,'XScale','log'); grid on;
    xlabel('centre frequency (Hz)'); ylabel('filter #');
    title('VoicesPerOctave: DENSITY'); legend(compose('%d', voices), 'Location','best');

    subplot(1,3,2); hold on;
    wavelist = {'Morse','amor','bump'};
    for i = 1:3
        fb = cwtfilterbank('SignalLength',N,'SamplingFrequency',fs,'Wavelet',wavelist{i});
        H  = freqz(fb);
        fr = linspace(0, fs/2, size(H,2));
        plot(fr, H(round(size(H,1)/2), :), 'LineWidth', 1.1);
    end
    grid on; xlim([0 300]); xlabel('frequency (Hz)'); ylabel('gain');
    title('Wavelet: SHAPE (middle filter)'); legend(wavelist, 'Location','best');

    subplot(1,3,3); hold on;
    lims = {[20 400], [50 150]};
    for i = 1:2
        fb = cwtfilterbank('SignalLength',N,'SamplingFrequency',fs,'FrequencyLimits',lims{i});
        plot(centerFrequencies(fb), 1:numel(scales(fb)), '.-', 'MarkerSize', 8);
    end
    set(gca,'XScale','log'); grid on;
    xlabel('centre frequency (Hz)'); ylabel('filter #');
    title('FrequencyLimits: REACH'); legend({'20-400 Hz','50-150 Hz'}, 'Location','best');
    sgtitle('cwtfilterbank -- Fig 1: design-stage parameters');
    fbd_save(f1, outDir, 'cwt_1_design');

    % ---------- Fig 2: filter properties ----------
    fb  = cwtfilterbank('SignalLength',N,'SamplingFrequency',fs);
    cf  = centerFrequencies(fb);  sc = scales(fb);
    pbw = powerbw(fb);            q  = qfactor(fb);

    f2 = fbd_fig(doExport, [70 70 1250 380]);
    subplot(1,3,1); semilogy(sc, '.-'); grid on;
    xlabel('filter #'); ylabel('scale'); title('scales (log-spaced)');
    subplot(1,3,2); semilogy(cf, '.-'); grid on;
    xlabel('filter #'); ylabel('centre frequency (Hz)'); title('centerFrequencies (Hz)');
    subplot(1,3,3);
    loglog(pbw.Frequencies, pbw.HalfPowerBandwidth, '.-'); grid on;
    xlabel('centre frequency (Hz)'); ylabel('half-power bandwidth (Hz)');
    title(sprintf('powerbw  |  qfactor = %.3f (ONE number)', q));
    sgtitle('cwtfilterbank -- Fig 2: filter properties. A CWT is CONSTANT-Q by construction');
    fbd_save(f2, outDir, 'cwt_2_properties');

    % ---------- Fig 3: signals with known features ----------
    fLo = 20;  fHi = 200;
    finst  = fLo + (fHi - fLo) * t;                       % instantaneous frequency
    xchirp = cos(2*pi*(fLo*t + 0.5*(fHi-fLo)*t.^2));
    fBurst = 100;  tBurst = 0.5;  sBurst = 0.02;
    xburst = exp(-(t - tBurst).^2 / (2*sBurst^2)) .* cos(2*pi*fBurst*t);
    x = xchirp + xburst;

    f3 = fbd_fig(doExport, [80 80 1150 420]);
    subplot(3,1,1); plot(t, xchirp); ylabel('chirp'); grid on;
    title(sprintf('linear chirp %g \\rightarrow %g Hz', fLo, fHi));
    subplot(3,1,2); plot(t, xburst); ylabel('burst'); grid on;
    title(sprintf('Gaussian burst: %g Hz at t = %g s', fBurst, tBurst));
    subplot(3,1,3); plot(t, x); ylabel('sum'); xlabel('time (s)'); grid on;
    sgtitle('cwtfilterbank -- Fig 3: signals with known features');
    fbd_save(f3, outDir, 'cwt_3_signal');

    % ---------- Fig 4: decomposition ----------
    [cfs, fcwt] = wt(fb, x);
    A = abs(cfs);
    yt = round(linspace(1, numel(fcwt), 6));

    f4 = fbd_fig(doExport, [90 90 1250 720]);
    subplot(2,2,[1 2]);
    imagesc(t, 1:numel(fcwt), A); set(gca,'YDir','normal'); hold on;
    [~, iTrue] = min(abs(fcwt(:) - finst), [], 1);
    plot(t, iTrue, 'w--', 'LineWidth', 1.2);
    set(gca,'YTick',yt,'YTickLabel',compose('%.0f', fcwt(yt)));
    ylabel('frequency (Hz)'); xlabel('time (s)');
    title('scalogram |wt|  (dashed = true instantaneous frequency)');
    subplot(2,2,3); semilogx(fcwt, timeSpectrum(fb, x), '.-'); grid on;
    xlabel('frequency (Hz)'); title('timeSpectrum (time-averaged \rightarrow per scale)');
    subplot(2,2,4); plot(t, scaleSpectrum(fb, x)); grid on;
    xlabel('time (s)'); title('scaleSpectrum (scale-averaged \rightarrow per time)');
    sgtitle('cwtfilterbank -- Fig 4: decomposition and marginals');
    fbd_save(f4, outDir, 'cwt_4_decomposition');

    % ---------- Fig 5: quantitative recovery ----------
    % ⚠ Measured over the INTERIOR only. The wavelet's time support makes the record edges
    % unreliable -- the cone of influence -- so the guard is stated rather than silently
    % trimmed. This is the CWT's counterpart of graphfilterbank's MassLost.
    guard = round(0.15*N);
    ii = guard:(N-guard);
    [~, ipk] = max(A(:, ii), [], 1);
    fRidge = fcwt(ipk).';
    ridgeErr = median(abs(fRidge - finst(ii)));

    % Restricted to a window around the burst, so the chirp's own crossing of 100 Hz
    % cannot be mistaken for it.
    win = abs(t - tBurst) < 4*sBurst;
    Aw = A;  Aw(:, ~win) = 0;
    [~, imax] = max(Aw(:));
    [iF, iT] = ind2sub(size(Aw), imax);

    % Tolerances come from the bank's own resolution, not from taste.
    tolF = interp1(pbw.Frequencies, pbw.HalfPowerBandwidth, fBurst, 'linear', 'extrap');
    checks = [ ...
        fbd_check('chirp ridge median error', ridgeErr, 0, 0.06*fHi, 'Hz'), ...
        fbd_check('burst centre frequency',   fcwt(iF), fBurst, tolF, 'Hz'), ...
        fbd_check('burst centre time',        t(iT), tBurst, 4*sBurst, 's')];

    f5 = fbd_fig(doExport, [100 100 1100 420]);
    subplot(1,2,1);
    plot(t(ii), finst(ii), 'k-', 'LineWidth', 1.3); hold on;
    plot(t(ii), fRidge, 'r.', 'MarkerSize', 5); grid on;
    xlabel('time (s)'); ylabel('frequency (Hz)');
    legend({'true','measured ridge'}, 'Location','northwest');
    title(sprintf('chirp: median error %.2f Hz', ridgeErr));
    subplot(1,2,2);
    imagesc(t, 1:numel(fcwt), Aw); set(gca,'YDir','normal'); hold on;
    plot(t(iT), iF, 'ro', 'MarkerSize', 11, 'LineWidth', 1.6);
    set(gca,'YTick',yt,'YTickLabel',compose('%.0f', fcwt(yt)));
    xlabel('time (s)'); ylabel('frequency (Hz)');
    title(sprintf('burst found at %.1f Hz, %.3f s', fcwt(iF), t(iT)));
    sgtitle('cwtfilterbank -- Fig 5: measured vs analytic');
    fbd_save(f5, outDir, 'cwt_5_recovery');

    [ok, ~] = fbd_report(checks);
    out.checks = checks;  out.ok = ok;
    out.ridgeErrHz = ridgeErr;  out.burstFreqHz = fcwt(iF);  out.burstTimeS = t(iT);
    out.qfactor = q;  out.numFilters = numel(sc);
end

% Author: Diellor Basha, 2026
