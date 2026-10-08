classdef tScaleCleanspan < matlab.unittest.TestCase
% TSCALECLEANSPAN  rheome.scale.cleanspan and the power-ratio normalisation periodicflow screens tiles with.
%
% The planted case is a real recording's: a ~3.5 s common-mode artefact ending 0.5 s before the last
% sample, and a hot first 0.5 s. No cached data.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function cleanRecordingIsNotTrimmed(tc)
            rng(1);  fs = 600;  F = randn(40, 60*fs);
            C = rheome.scale.cleanspan(F, fs);
            tc.verifyEqual([C.first C.last], [1 size(F,2)]);
            tc.verifyFalse(any(C.mask));
        end

        function endArtefactAndHotStartAreTrimmed(tc)
            rng(2);  fs = 600;  N = 127*fs;  F = randn(40, N);
            t = (123.5*fs + 1) : (126.5*fs);                      % that recording: 123.5-126.5 s, last 0.5 s clean
            F(:, t) = F(:, t) + 100 * randn(1, numel(t));          % common mode
            F(:, 1:fs/2) = 10 * F(:, 1:fs/2);
            C = rheome.scale.cleanspan(F, fs);
            tc.verifyGreaterThan(C.first, fs/2);                   % past the hot start
            tc.verifyLessThanOrEqual(C.first, 1.5*fs + 1);         % ... by at most the dilation
            tc.verifyLessThan(C.last, 123.5*fs + 1);               % before the artefact
            tc.verifyGreaterThanOrEqual(C.last, 122.5*fs);
            tc.verifyFalse(any(C.mask(C.first:C.last)));
        end

        function interiorArtefactIsMaskedNotTrimmed(tc)
            rng(3);  fs = 600;  F = randn(40, 120*fs);
            t = 60*fs + (1:fs);  F(:, t) = 50 * F(:, t);
            C = rheome.scale.cleanspan(F, fs);
            tc.verifyEqual([C.first C.last], [1 size(F,2)]);
            tc.verifyTrue(all(C.mask(t)));
            tc.verifyFalse(any(C.mask([1:58*fs, 63*fs:end])));
        end

        function slowOnsetIsCaughtByTheLowBand(tc)
            % its tail turns red before it turns loud: a 0.3 Hz common-mode swell at ~3x the
            % low-band level that barely moves the broadband std must still be trimmed.
            rng(5);  fs = 600;  N = 120*fs;  F = randn(40, N);
            t = (110*fs + 1) : N;  tt = (0:numel(t)-1) / fs;
            [b, a] = butter(4, 4/(fs/2), 'low');  lowStd = median(std(filtfilt(b, a, F.'), 0, 1));
            F(:, t) = F(:, t) + 6 * lowStd * sin(2*pi*0.3*tt);
            C0 = rheome.scale.cleanspan(F, fs, LowHz=0);
            C  = rheome.scale.cleanspan(F, fs);
            tc.verifyEqual(C0.last, N);                            % broadband alone misses it
            tc.verifyLessThanOrEqual(C.last, 110*fs);
        end

        function rectangularFftOverWelchIsFsNtOverTwo(tc)
            % periodicflow's power_ratio divides the data-derived scale by fs*nT/2: for a one-sided
            % Welch PSD and an unwindowed FFT that is the expected |X_k|^2 / P(f_k).
            rng(4);  fs = 600;  F = randn(30, 120*fs);
            nf = 2^nextpow2(8*fs);  [P, f] = pwelch(F.', hann(nf), nf/2, nf, fs);
            nT = 6*fs;  W = fft(F(:, 1:nT), [], 2);  fw = (0:nT-1)*fs/nT;  pos = 2:nT/2;
            Pi = interp1(f, P, fw(pos)')';  in = fw(pos) >= 1 & fw(pos) <= 45;
            sc = mean(abs(W(:, pos(in))).^2, 'all') / mean(Pi(:, in), 'all');
            tc.verifyEqual(sc / (fs*nT/2), 1, 'AbsTol', 0.1);
        end
    end
end

% Author: Diellor Basha, 2026
