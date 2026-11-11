classdef tTimeFilterBank < matlab.unittest.TestCase
% The designed constant-Q frame: tight on the whole axis by construction, compactly
% supported, evaluated per member at its own rate with exact magnitudes, Parseval-exact,
% and anchored at an absolute frequency so every rate shares the grid.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        N  = 4000
        tfb
    end

    methods (TestClassSetup)
        function make(tc)
            tc.tfb = rheome.timefilterbank(tc.N, 'SamplingFrequency', tc.fs);
        end
    end

    methods (Test)

        function theFrameIsTightOnTheWholeAxis(tc)
            b = framebounds(tc.tfb);
            tc.verifyEqual(b.A, 2, 'AbsTol', 1e-12);                    % Gain^2/2 for the default Gain = 2
            tc.verifyEqual(b.B, 2, 'AbsTol', 1e-12);
            tc.verifyEqual(b.Uncovered, 0);
            tc.verifyTrue(isframetight(tc.tfb));
        end

        function everyMemberIsZeroOutsideItsBins(tc)
            [H, ~] = freqz(tc.tfb);  B = bins(tc.tfb);
            for m = 1:tc.tfb.NumMembers
                out = true(1, size(H, 2));  out(B(m,1):B(m,2)) = false;
                tc.verifyEqual(max(abs(H(m, out))), 0, sprintf('member %d', m));
                tc.verifyEqual(max(H(m, :)), tc.tfb.Gain, 'RelTol', 2e-3);   % the sampled peak IS the gain (2): the grid misses the centre by < 0.1%
            end
            k = kinds(tc.tfb);
            tc.verifyEqual(k{1}, 'lowpass');  tc.verifyEqual(k{end}, 'highpass');
            tc.verifyTrue(all(strcmp(k(2:end-1), 'band')));
        end

        function membersSitOnTheAnchoredGridAndFormOctaves(tc)
            fc = centerFrequencies(tc.tfb);
            V = tc.tfb.VoicesPerOctave;
            k = round(V * log2(fc(2:end-1)));                         % anchor 1 Hz
            tc.verifyEqual(fc(2:end-1), 2.^(k / V), 'RelTol', 1e-12);
            tc.verifyEqual(diff(k), ones(1, numel(k)-1));              % consecutive voices
        end

        function subBandMagnitudesAreExactAtCoincidentInstants(tc)
            % Force Nm = N for a mid member by a bank whose bins fit: compare |coef| with
            % |ifft(fft(x).*H)| where the instants coincide (N/Nm integer).
            rng(3);  x = randn(tc.N, 1);
            [H, ~] = freqz(tc.tfb);  X = fft(x);
            C = wt(tc.tfb, x);
            for m = [2, round(tc.tfb.NumMembers/2), tc.tfb.NumMembers-1]
                full = ifft(X .* [H(m, :), zeros(1, tc.N - size(H, 2))].');
                Nm = numel(C(m).coef);  r = tc.N / Nm;
                if r == round(r)                                       % coincident instants
                    idx = 1:r:tc.N;
                    tc.verifyEqual(abs(C(m).coef), abs(full(idx)), 'RelTol', 1e-10, 'AbsTol', 1e-12);
                else                                                   % interpolated: compare energy
                    tc.verifyEqual(sum(abs(C(m).coef).^2) * C(m).weight, sum(abs(full).^2), 'RelTol', 1e-10);
                end
            end
        end

        function parsevalHoldsAcrossMembers(tc)
            rng(4);  x = randn(tc.N, 1);
            C = wt(tc.tfb, x);
            E = sum(arrayfun(@(c) sum(abs(c.coef).^2) * c.weight, C));
            X = fft(x);
            excluded = (abs(X(1))^2 + abs(X(tc.N/2+1))^2) / tc.N;      % DC and Nyquist bins
            tc.verifyEqual(E, tc.tfb.Gain^2/2 * (sum(x.^2) - excluded), 'RelTol', 1e-10);   % (Gain^2/2) x the record's energy
        end

        function oversampleChangesLengthsNotEnergies(tc)
            rng(5);  x = randn(tc.N, 1);
            t2 = rheome.timefilterbank(tc.N, 'SamplingFrequency', tc.fs, 'Oversample', 2);
            t4 = rheome.timefilterbank(tc.N, 'SamplingFrequency', tc.fs, 'Oversample', 4);
            C2 = wt(t2, x);  C4 = wt(t4, x);
            for m = 1:t2.NumMembers
                tc.verifyGreaterThanOrEqual(numel(C4(m).coef), numel(C2(m).coef));
                tc.verifyEqual(sum(abs(C4(m).coef).^2) * C4(m).weight, sum(abs(C2(m).coef).^2) * C2(m).weight, 'RelTol', 1e-10);
            end
        end

        function twoRatesShareTheGridExactly(tc)
            a = rheome.timefilterbank(6000,  'SamplingFrequency', 600);
            b = rheome.timefilterbank(24000, 'SamplingFrequency', 2400);
            fa = centerFrequencies(a);  fb = centerFrequencies(b);
            fa = fa(2:end-1);  fb = fb(2:end-1);
            common = intersect(round(1e9*fa), round(1e9*fb)) / 1e9;
            tc.verifyGreaterThan(numel(common), 20);
            tc.verifyEqual(numel(common), numel(fa));                   % every 600 Hz member is a 2400 Hz member
        end

        function explicitLimitsBoundTheBandMembersOnly(tc)
            t = rheome.timefilterbank(tc.N, 'SamplingFrequency', tc.fs, 'FrequencyLimits', [2 40]);
            fc = centerFrequencies(t);  fc = fc(2:end-1);
            tc.verifyGreaterThanOrEqual(min(fc), 2);
            tc.verifyLessThanOrEqual(max(fc) * 2^(1/t.VoicesPerOctave), 40 * (1 + 1e-9));
            b = framebounds(t);
            tc.verifyEqual(b.A, 2, 'AbsTol', 1e-12);                   % edge members still complete it
        end

        function supportIsScaleCovariantAndQIsConstant(tc)
            s = support(tc.tfb);
            tc.verifyGreaterThan(s.timeTimesFc, 0.5);
            tc.verifyLessThan(s.timeTimesFc, 50);
            ok = s.perMember(1, :) > 0;
            tc.verifyLessThan(std(s.perMember(1, ok)) / median(s.perMember(1, ok)), 0.05);
            tc.verifyGreaterThan(s.Q, 3);  tc.verifyLessThan(s.Q, 20);
        end

        function gainIsThePeakAndAToneOfAmplitudeAGivesGHalfA(tc)
            % one meaning of gain: the peak in-band response, so |c| = (Gain/2) * A for a real tone on a member's centre
            A = 1.5;  t = (0:tc.N-1)' / tc.fs;
            for G = [1 2 3]
                tf = rheome.timefilterbank(tc.N, 'SamplingFrequency', tc.fs, 'Gain', G);
                tc.verifyEqual(tf.Gain, G);
                [H, ~] = freqz(tf);
                tc.verifyEqual(max(H, [], 2), G * ones(tf.NumMembers, 1), 'RelTol', 2e-3, sprintf('sampled peak, Gain %g', G));
                fcb = centerFrequencies(tf);                            % a band member's exact peak is AT its centre
                tc.verifyEqual(gains(tf, 5, fcb(5)), G, 'AbsTol', 1e-12);
                b = framebounds(tf);                                    % tight at any gain: Gain^2/2
                tc.verifyEqual([b.A b.B], G^2/2 * [1 1], 'AbsTol', 1e-12);
                C = wt(tf, A * cos(2*pi*8*t));                          % 8 Hz = 2^(12/4) is a member centre, on a bin
                fc = centerFrequencies(tf);  [~, m] = min(abs(fc - 8));
                tc.verifyEqual(max(abs(C(m).coef)), G/2 * A, 'RelTol', 1e-9, sprintf('|c|, Gain %g', G));
            end
            tc.verifyError(@() rheome.timefilterbank(tc.N, 'SamplingFrequency', tc.fs, 'Gain', 0), ?MException);
        end

        function aToneLandsInTheMembersAroundIt(tc)
            f0 = 10;  t = (0:tc.N-1)' / tc.fs;  x = sin(2*pi*f0*t);
            C = wt(tc.tfb, x);  fc = centerFrequencies(tc.tfb);
            E = arrayfun(@(c) sum(abs(c.coef).^2) * c.weight, C);
            near = abs(log2(fc / f0)) <= 2 / tc.tfb.VoicesPerOctave;
            tc.verifyGreaterThan(sum(E(near)) / sum(E), 0.999);
        end

    end
end
