classdef tSpectralSplit < matlab.unittest.TestCase
% rheome.spectral.gains and rheome.spectral.decompose, the periodic / aperiodic split every MS1 split result
% rests on (G4, G10): the power partition is exact, clipping keeps it exact, a 1/f background without a rhythm
% leaves (almost) nothing periodic, a planted rhythm is kept, and the unit trap is caught.
%
% Author: Diellor Basha, 2026

    properties (Constant)
        fs = 200
        nT = 200 * 6           % one 6 s window, as the G10 split
    end

    methods (Static)
        function [Fw, f, Pfit, ff] = i_window(X, fs)
            % the window's positive-frequency FFT and a Welch PSD in |FFT bin|^2 units (decompose's contract)
            nT = size(X, 2);  Fw = fft(X, [], 2);  f = (0:nT-1)*fs/nT;  pos = 2:floor(nT/2)+1;
            [P, ff] = pwelch(X.', hann(2*fs), fs, 2*fs, fs);  m = ff >= 1 & ff <= 45;
            Fw = Fw(:, pos);  f = f(pos)';  Pfit = (P(m, :) * fs*nT/2)';  ff = ff(m);
        end
    end

    methods (Test)

        function gainsArePowerComplementary(tc)
            rng(1);  P = rand(50, 3) + 0.1;  Pap = P .* rand(50, 3);
            [hp, ha] = rheome.spectral.gains(P, Pap);
            tc.verifyEqual(hp.^2 + ha.^2, ones(50, 3), 'AbsTol', 1e-14);
            tc.verifyEqual(ha, sqrt(Pap ./ P), 'AbsTol', 1e-14);
            tc.verifyTrue(all(hp >= 0 & hp <= 1 & ha >= 0 & ha <= 1, 'all'));
        end

        function anOvershootIsClippedAndStillPartitions(tc)
            % the fit above the data: all aperiodic there, and the identity holds
            P = [1; 2; 3];  Pap = [2; 1; 3];
            [hp, ha] = rheome.spectral.gains(P, Pap);
            tc.verifyEqual([hp ha], [0 1; sqrt(0.5) sqrt(0.5); 0 1], 'AbsTol', 1e-14);
            [hp0, ha0] = rheome.spectral.gains(P, zeros(3, 1));
            tc.verifyEqual([hp0 ha0], [ones(3, 1) zeros(3, 1)]);
        end

        function gainsRejectMismatchedSizes(tc)
            tc.verifyError(@() rheome.spectral.gains(ones(4, 2), ones(2, 4)), 'spectral:gains:size');
        end

        function decomposePartitionsPowerExactly(tc)
            rng(2);  X = rheome.spectral.synthesise(tc.nT, tc.fs, 1.5, struct('nS', 4));
            [Fw, f, Pfit, ff] = tc.i_window(X, tc.fs);
            D = rheome.spectral.decompose(Fw, f, struct('fitPower', Pfit, 'fitF', ff, 'knee', true));
            tc.verifyLessThan(D.partition, 1e-12);
            tc.verifyEqual(D.Cper + D.Cap, (D.hper + D.hap) .* Fw, 'AbsTol', 1e-9 * max(abs(Fw(:))));
            tc.verifySize(D.hper, size(Fw));
            tc.verifyEqual(D.ap.exponent, 1.5*ones(1, 4), 'AbsTol', 0.4);
        end

        function aPlantedRhythmIsKeptAndTheBackgroundIsNot(tc)
            % 1/f + a 10 Hz sinusoid: the 10 Hz bin is periodic, the 30-40 Hz background is not
            rng(3);  X = rheome.spectral.synthesise(tc.nT, tc.fs, 1.5, struct('nS', 3));
            t = (0:tc.nT-1)/tc.fs;  X = X + 3*std(X(:))*sin(2*pi*10*t);
            [Fw, f, Pfit, ff] = tc.i_window(X, tc.fs);
            D = rheome.spectral.decompose(Fw, f, struct('fitPower', Pfit, 'fitF', ff, 'knee', true));
            share = @(b) sum(abs(D.Cper(:, b)).^2, 'all') / sum(abs(Fw(:, b)).^2, 'all');
            tc.verifyGreaterThan(share(abs(f - 10) < 0.01), 0.95);
            tc.verifyLessThan(share(f >= 30 & f <= 40), 0.5);
        end

        function decomposeCatchesTheUnitTrap(tc)
            % a PSD passed without the |FFT bin|^2 conversion sits far above the coefficients
            rng(4);  X = rheome.spectral.synthesise(tc.nT, tc.fs, 1.5, struct('nS', 2));
            [Fw, f, Pfit, ff] = tc.i_window(X, tc.fs);
            tc.verifyError(@() rheome.spectral.decompose(Fw, f, struct('fitPower', Pfit * 1e6, 'fitF', ff)), ...
                           'spectral:decompose:scale');
            tc.verifyError(@() rheome.spectral.decompose(Fw, f, struct('fitPower', Pfit)), 'spectral:decompose:fitF');
            tc.verifyError(@() rheome.spectral.decompose(Fw, f, struct('fitPower', Pfit(1, :), 'fitF', ff)), ...
                           'spectral:decompose:fitRows');
            tc.verifyError(@() rheome.spectral.decompose(Fw, f(1:end-1)), 'spectral:decompose:size');
        end
    end
end

% Author: Diellor Basha, 2026
