function sp = split(Y, fs, mask, opts)
% MODES.SPLIT  Periodic and aperiodic parts of tile series, by the Welch background fit of each row.
%
%   sp = rheome.modes.split(Y, fs)                 % Y [nTiles x N] real, fs Hz
%   sp = rheome.modes.split(Y, fs, mask, FitRange=[1 45], SegS=8)
%
% The background of each row is fitted with a knee on its Welch PSD (8 s Hann, 50 %, segments touching an
% artefact block in mask left out: rheome.scale.cleanwelch; rheome.spectral.aperiodic), and the complementary
% gains of rheome.spectral.gains, hper = sqrt(max(0, 1 - Pap/P)) and hap = sqrt(Pap/P), are applied to the
% row's whole-span FFT as zero-phase filters. Outside FitRange nothing is called rhythm: hper = 0, hap = 1.
% The gains come from the Welch PSD, not the span's own periodogram, which is biased low in log power
% (rheome.spectral.decompose, "fit on a Welch PSD").
%
%   sp.per sp.ap  [nTiles x N]  the periodic and aperiodic series (power partitioned in expectation)
%   sp.f [nF x 1]  sp.P sp.Pap [nF x nTiles]  Welch PSD and fitted background on the fit range
%   sp.exponent [1 x nTiles]
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        fs (1,1) double {mustBePositive}
        mask = []
        opts.FitRange (1,2) double = [1 45]
        opts.SegS (1,1) double = 8
    end
    N = size(Y, 2);  if isempty(mask), mask = false(1, N); end
    nf = 2^nextpow2(opts.SegS * fs);
    if N < nf, error('modes:split:short', 'span of %.1f s is shorter than one %g s Welch segment', N/fs, nf/fs); end
    [P, f] = rheome.scale.cleanwelch(Y, fs, 1:N, logical(mask(:)'), nf);
    mf = f >= opts.FitRange(1) & f <= opts.FitRange(2);
    ap = rheome.spectral.aperiodic(P(mf, :), f(mf), struct('knee', true));
    [hp, ha] = rheome.spectral.gains(P(mf, :), ap.Pap);
    fk = (0:N-1)' * fs / N;  pos = 2:floor(N/2) + 1;
    Hp = zeros(N, size(Y, 1));  Ha = ones(N, size(Y, 1));
    Hp(pos, :) = interp1(f(mf), hp, fk(pos), 'linear', 0);
    Ha(pos, :) = interp1(f(mf), ha, fk(pos), 'linear', 1);
    neg = N - pos + 2;  ok = neg > pos(end);                     % mirror onto the negative frequencies
    Hp(neg(ok), :) = Hp(pos(ok), :);  Ha(neg(ok), :) = Ha(pos(ok), :);
    Fy = fft(Y, [], 2);
    sp.per = real(ifft(Fy .* Hp.', [], 2));  sp.ap = real(ifft(Fy .* Ha.', [], 2));
    sp.f = f(mf);  sp.P = P(mf, :);  sp.Pap = ap.Pap;  sp.exponent = ap.exponent;
end

% Author: Diellor Basha, 2026
