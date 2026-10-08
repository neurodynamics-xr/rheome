classdef tSpectralSynth < matlab.unittest.TestCase
    % Round trip between rheome.spectral.synthesise and rheome.spectral.aperiodic.
    %
    % Author: Diellor Basha, 2026

    properties (Constant)
        fs = 200
        nT = 200 * 60          % 60 s, the length the background synthesis actually uses
    end

    methods (Static)
        function [chiHat, r2] = i_fit(x, fs, knee)
            [P,f] = pwelch(x(:), hann(round(4*fs)), [], [], fs);
            keep  = f > 1 & f < 45;
            ap    = rheome.spectral.aperiodic(P(keep), f(keep), struct('knee', knee));
            chiHat = ap.exponent;  r2 = ap.r2;
        end
    end

    methods (Test)

        function theSeriesIsReal(tc)
            % ⚠ the whole point of the Hermitian mirror
            rng(1);
            X = rheome.spectral.synthesise(tc.nT, tc.fs, 2.0);
            tc.verifyTrue(isreal(X), 'synthesise returned a complex series.');
            tc.verifySize(X, [1 tc.nT]);
        end

        function anOddLengthIsAlsoReal(tc)
            % no Nyquist bin to special-case; the mirror indices must still line up
            rng(2);
            X = rheome.spectral.synthesise(1001, tc.fs, 1.5);
            tc.verifyTrue(isreal(X));
            tc.verifySize(X, [1 1001]);
        end

        function plantedExponentIsRecovered(tc)
            % the round trip that licenses using this as a background
            rng(3);
            for chi = [0.5 1.0 2.0 3.0]
                X = rheome.spectral.synthesise(tc.nT, tc.fs, chi);
                [chiHat, r2] = tc.i_fit(X, tc.fs, false);
                tc.verifyEqual(chiHat, chi, 'AbsTol', 0.15, ...
                    sprintf('planted chi %.2f came back %.2f.', chi, chiHat));
                % ⚠ a shallow spectrum scores lower: chi 0.5 fits at r2 ~ 0.88 because the
                %   scatter is a larger share of a smaller dynamic range
                tc.verifyGreaterThan(r2, 0.85);
            end
        end

        function aKneeFlattensTheLowEnd(tc)
            % a knee must lower the apparent no-knee exponent, since it flattens below k^(1/chi)
            rng(4);
            plain = rheome.spectral.synthesise(tc.nT, tc.fs, 2.5);
            kneed = rheome.spectral.synthesise(tc.nT, tc.fs, 2.5, struct('knee', 300));
            tc.verifyLessThan(tc.i_fit(kneed, tc.fs, false), tc.i_fit(plain, tc.fs, false));
        end

        function rowsAreIndependentAndUnitVariance(tc)
            rng(5);
            X = rheome.spectral.synthesise(tc.nT, tc.fs, 1.8, struct('nS', 6));
            tc.verifySize(X, [6 tc.nT]);
            tc.verifyEqual(std(X, 0, 2), ones(6,1), 'AbsTol', 1e-10);
            % ⚠ correlate the DIFFERENCED series. Independent 1/f draws correlate at up to
            %   0.65 in the raw series -- they are dominated by a few low-frequency bins, so
            %   the effective DOF is tiny. Differencing whitens, and the same rows give 0.03.
            D = abs(corr(diff(X,1,2).'));  D(1:7:end) = 0;
            tc.verifyLessThan(max(D(:)), 0.15, 'rows are not independent.');
        end

        function inBandPowerIsMatchedExactly(tc)
            % the level knob, which is what lets a background sit under planted activity.
            % ⭐ matched on the REALISATION, so this is exact, not statistical
            rng(6);
            want = [4; 25];
            X = rheome.spectral.synthesise(tc.nT, tc.fs, 2.0, ...
                struct('nS', 2, 'power', want, 'band', [1 45]));
            got = zeros(2,1);
            for i = 1:2
                W = abs(fft(X(i,:))).^2 / tc.nT^2;
                k = 1:floor(tc.nT/2);  f = k * tc.fs / tc.nT;
                got(i) = 2 * sum(W(k(f >= 1 & f <= 45)));
            end
            tc.verifyEqual(got, want, 'RelTol', 1e-10, ...
                sprintf('in-band power %.4f/%.4f against %.4f/%.4f.', got(1), got(2), want(1), want(2)));
        end

        function aWelchEstimateSeesRoughlyTheRequestedPower(tc)
            % the caller measures power with pwelch, so the Parseval definition must not be
            % far from what a Welch fit sees. ⚠ they do not agree exactly: a 4 s hann window
            % leaks the large sub-1 Hz content into the bottom of the band.
            rng(8);
            want = 10;
            x = rheome.spectral.synthesise(tc.nT, tc.fs, 2.0, struct('power', want, 'band', [1 45]));
            [P,f] = pwelch(x(:), hann(round(4*tc.fs)), [], [], tc.fs);
            got = sum(P(f > 1 & f < 45)) * mean(diff(f));
            tc.verifyEqual(got / want, 1, 'RelTol', 0.5, ...
                sprintf('welch in-band power %.3f against requested %.3f.', got, want));
        end

        function shapeIsIndependentOfLevel(tc)
            % ⭐ chi and knee fix the shape; opts.power fixes the level
            rng(7);
            a = rheome.spectral.synthesise(tc.nT, tc.fs, 2.0);
            rng(7);
            b = rheome.spectral.synthesise(tc.nT, tc.fs, 2.0, struct('power', 1e6, 'band', [1 45]));
            tc.verifyEqual(tc.i_fit(b, tc.fs, false), tc.i_fit(a, tc.fs, false), 'AbsTol', 1e-6);
        end

        function aMeasuredShapeIsReproduced(tc)
            % ⭐ the path that exists because the knee grid tops out at 1000 and real MEG pins it
            f  = (1:0.25:45)';
            Pt = 1./(1000 + f.^2.7);
            X  = rheome.spectral.synthesise(tc.nT, tc.fs, NaN, ...
                     struct('nS', 40, 'shape', Pt, 'shapeF', f));
            [P, ff] = pwelch(X.', hann(round(4*tc.fs)), [], [], tc.fs);
            P = mean(P, 2);  k = ff > 1 & ff < 45;
            Pi = 10.^interp1(log10(f), log10(Pt), log10(ff(k)), 'linear', 'extrap');
            tc.verifyGreaterThan(corr(log(P(k)), log(Pi)), 0.99);
            % ⭐ band-to-band faithfulness: the ratio must be CONSTANT, level being arbitrary
            r = zeros(1,4);  bnd = [2 4; 4 8; 8 16; 16 32];
            for i = 1:4
                kb = ff >= bnd(i,1) & ff <= bnd(i,2) & k;
                ki = ff(k) >= bnd(i,1) & ff(k) <= bnd(i,2);
                r(i) = sum(P(kb))/sum(Pi(ki));
            end
            tc.verifyLessThan(max(r)/min(r), 1.1, 'the shape is not faithful band to band');
        end

        function aShapeNeedsItsFrequencies(tc)
            tc.verifyError(@() rheome.spectral.synthesise(tc.nT, tc.fs, NaN, ...
                struct('shape', [1;2;3])), ?MException);
            tc.verifyError(@() rheome.spectral.synthesise(tc.nT, tc.fs, NaN, ...
                struct('shape', [1;2;3], 'shapeF', [1;2])), ?MException);
        end

        function badInputsAreRejected(tc)
            tc.verifyError(@() rheome.spectral.synthesise(4, tc.fs, 2.0), ?MException);
            tc.verifyError(@() rheome.spectral.synthesise(tc.nT, tc.fs, NaN), ?MException);
            tc.verifyError(@() rheome.spectral.synthesise(tc.nT, tc.fs, 2.0, ...
                struct('nS', 3, 'power', [1;2])), ?MException);
        end

    end
end

% Author: Diellor Basha, 2026
