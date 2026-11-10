classdef tVortexSignal < matlab.unittest.TestCase
    % rheome.cores burst stages, band frames and the lag-free envelope-matched null against nxr's synthetic goldens
    % (tests/goldens/vortex-signal, nxr-analysis commit 30b9622). nxr is the reference.
    %
    % Float32 stages are compared at 1e-5 of the peak. The null cannot reuse numpy's generator: it is fed the stored draw
    % (null_white). The risefall fits use lsqcurvefit against scipy's trf, so they agree to a tolerance, not exactly.
    %
    % Author: Diellor Basha, 2026

    properties
        G
        P
        Fit   % rheome.cores.burststages result
    end

    methods (TestClassSetup)
        function load(t)
            d = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'vortex-signal');
            M = jsondecode(fileread(fullfile(d, 'manifest.json')));
            nm = fieldnames(M.arrays);
            for i = 1:numel(nm), t.G.(nm{i}) = tVortexSignal.readArray(d, M.arrays.(nm{i})); end
            t.P = M.parameters;
        end
    end

    methods (Static)
        function x = readArray(d, a)
            sh = double(a.shape(:)).';
            f = fopen(fullfile(d, a.file), 'r'); x = fread(f, prod(sh), ['*' strrep(a.dtype, 'float64', 'double')]); fclose(f);
            if numel(sh) == 1, x = x(:); else, x = permute(reshape(x, fliplr(sh)), numel(sh):-1:1); end
        end

        function e = relpeak(x, ref)
            e = max(abs(double(x(:)) - ref(:))) / max(abs(ref(:)));
        end

        function z = frames(t)
            z = complex(single(t.G.frames_real), single(t.G.frames_imag));
        end
    end

    methods (Test)
        function morletGain(t)
            f = t.G.fft_freq_hz(:).';
            g = rheome.cores.morletgain(f, t.P.f0_hz, t.P.fwhm_hz);
            t.verifyLessThan(max(abs(g - t.G.morlet_gain(:).')), 1e-14, 'morlet gain');
            t.verifyEqual(max(g), 2, AbsTol=1e-12);
            t.verifyEqual(g(f <= 0), zeros(1, nnz(f <= 0)));
        end

        function bandFrames(t)
            z = rheome.cores.bandframes(t.G.sensors, t.P.fs_hz, t.P.f0_hz, t.P.fwhm_hz, t.P.frames_per_s);
            t.verifyClass(z, 'single');
            t.verifySize(z, [6 5760]);
            pk = max(abs(t.G.frames_real(:)));
            t.verifyLessThan(tVortexSignal.relpeak(real(z), t.G.frames_real), 1e-5, 'frames real');
            t.verifyLessThan(tVortexSignal.relpeak(imag(z), t.G.frames_imag), 1e-5, 'frames imag');
            t.assertGreaterThan(pk, 1);
        end

        function bandedCosine(t)
            % band-limited A cos(2 pi f0 t + phi) -> A exp(i(...)): amplitude A, not A/2
            tt = (0:2999) / 300;  x = 3 * cos(2*pi*10*tt + 0.7);
            z = rheome.cores.bandframes(x, 300, 10, 4, 360);
            fr = (0:3599) / 360;
            t.verifyEqual(double(z), 3 * exp(1i * (2*pi*10*fr + 0.7)), AbsTol=1e-4);
        end

        function smoothing(t)
            x = [1 2 3 4 5 6 7 8];
            % even k: window [i-k/2, i+k/2-1], zero padded, divided by k
            y = rheome.cores.boxsmooth(x, 4);
            t.verifyEqual(y(1), (1 + 2) / 4, AbsTol=1e-14);          % x[-2..1] -> 1,2
            t.verifyEqual(y(4), (2 + 3 + 4 + 5) / 4, AbsTol=1e-14);   % x[1..4] 0-based -> 2..5
            t.verifyEqual(y(8), (6 + 7 + 8) / 4, AbsTol=1e-14);
            y3 = rheome.cores.boxsmooth(x, 3);
            t.verifyEqual(y3([1 5 8]), [3 15 15] / 3, AbsTol=1e-14);
        end

        function smoothedAmplitude(t)
            z = tVortexSignal.frames(t);
            a = sqrt(mean(abs(double(z)).^2, 1));
            t.verifyLessThan(tVortexSignal.relpeak(a, t.G.amp_mean), 1e-5, 'amp_mean');
            a1 = rheome.cores.boxsmooth(a, 36);
            t.verifyLessThan(tVortexSignal.relpeak(a1, t.G.amp_smooth1), 1e-5, 'amp_smooth1');
            ae = rheome.cores.boxsmooth(a1, 72);
            t.verifyLessThan(tVortexSignal.relpeak(ae, t.G.amp_smooth_env), 1e-5, 'amp_smooth_env');
        end

        function envelopeThresholdsAndCandidates(t)
            a1 = t.G.amp_smooth1(:).';  good = logical(t.G.good_frames(:).');
            [c, th] = rheome.cores.burstenvelope(a1, good, 36);
            t.verifyLessThan(tVortexSignal.relpeak(th, t.G.envelope_thresholds), 1e-12, 'thresholds');
            ref = double(t.G.candidates);
            t.verifyEqual(c, [ref(:,1) + 1, ref(:,2)], 'candidates (1-based inclusive)');
        end

        function percentileIsNumpyLinear(t)
            t.verifyEqual(rheome.cores.linpercentile([1 2 3 4], 50), 2.5);
            t.verifyEqual(rheome.cores.linpercentile([1 2 3 4 5], 5), 1.2, AbsTol=1e-14);
            t.verifyEqual(rheome.cores.linpercentile([3 1 2], 100), 3);
        end

        function oldBurstMask(t)
            z = tVortexSignal.frames(t);
            [k, info] = rheome.cores.burstmask(z, logical(t.G.good_frames(:).'), 36);
            t.verifyEqual(double(k), double(t.G.old_burst_mask(:).'), 'mask');
            t.verifyEqual(info.nBursts, 3);
            t.verifyEqual(info.burstCycles, 12.666666666666666, AbsTol=1e-12);
            t.verifyEqual(info.burstShare, 0.23870967741935484, AbsTol=1e-12);
        end

        function burstFitsAndMasks(t)
            z = tVortexSignal.frames(t);  good = logical(t.G.good_frames(:).');
            S = rheome.cores.burststages(z, good, 360, 36);
            t.verifyEqual(size(S.candidates, 1), 4);
            ref = t.G.risefall_fits;  riv = t.G.risefall_intervals_s;
            fn = {'tOn','tOff','tauR','tauF','b','h','r2'};
            tol = [1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-4];   % tOn/tOff s; tau relative; b,h relative; r2 absolute
            % achieved agreement is printed so the PR can state it
            worst = zeros(4, 7);
            for c = 1:4
                R = S.all(c);
                t.verifyTrue(R.ok); t.verifyEqual(double(R.atBound), ref(c, 9));
                for q = 1:7
                    d = abs(R.(fn{q}) - ref(c, q));
                    if q >= 3 && q <= 6, d = d / abs(ref(c, q)); end
                    worst(c, q) = d;
                    t.verifyLessThan(d, tol(q), sprintf('cand %d %s', c, fn{q}));
                end
                iv = [R.buildup R.plateau R.decay];
                t.verifyLessThan(max(abs(iv - riv(c, ~isnan(riv(c, :))))), 1e-3, sprintf('intervals %d', c));
            end
            disp(array2table(max(worst, [], 1), 'VariableNames', {'tOn_s','tOff_s','tauR_rel','tauF_rel','b_rel','h_rel','r2_abs'}));
            kept = t.G.burst_fits;
            t.verifySize(S.fits, size(kept));
            t.verifyLessThan(max(abs(S.fits(:,1:2) - kept(:,1:2)), [], 'all'), 1e-3, 'kept tOn/tOff');
            t.verifyLessThan(max(abs(S.fits(:,5) - kept(:,5))), 1e-4, 'kept r2');
            % masks: identical frame for frame at this tolerance (fits agree to ~1e-5 s, a frame is 2.8e-3 s)
            nm = {'onset', 'within', 'offset'};
            for q = 1:3
                ref = logical(t.G.(['mask_' nm{q}])(:).');
                t.verifyEqual(S.(nm{q}), ref, ['mask ' nm{q}]);
                t.verifyLessThan(nnz(S.(nm{q}) & ~good), 1);
            end
        end

        function nullStages(t)
            xr = single(t.G.frames_real);  good = logical(t.G.good_frames(:).');
            white = t.G.null_white;
            Y = rheome.cores.lagfreenull(xr, good, white);
            t.verifyClass(Y, 'single');
            t.verifyLessThan(tVortexSignal.relpeak(Y, t.G.null_lagfree), 1e-5, 'lagfree');
            E = rheome.cores.envelopenull(xr, good, 36, white);
            t.verifyClass(E, 'single');
            t.verifyLessThan(tVortexSignal.relpeak(E, t.G.null_envelope_matched), 1e-5, 'envelope matched');
        end

        function nullIntermediates(t)
            xr = double(single(t.G.frames_real));  good = logical(t.G.good_frames(:).');
            n = size(xr, 2);
            amp = sqrt(mean(abs(fft(xr, [], 2)).^2, 1));
            t.verifyLessThan(tVortexSignal.relpeak(amp(1:floor(n/2)+1), t.G.null_amp_spectrum), 1e-5, 'amp spectrum');
            A = cov(xr(:, good).');
            t.verifyLessThan(tVortexSignal.relpeak(A, t.G.null_cov), 1e-5, 'cov');
            [E, e] = eig((A + A.')/2, 'vector');  R = (E .* sqrt(max(e(:).', 0))) * E.';
            t.verifyLessThan(tVortexSignal.relpeak(R, t.G.null_sqrt_cov), 1e-5, 'sqrt cov');
            kk = 0:n-1;  full = amp(1:floor(n/2)+1);  full = full(min(kk, n - kk) + 1);
            W = ifft(fft(t.G.null_white, [], 2) .* full, [], 2, 'symmetric');
            t.verifyLessThan(tVortexSignal.relpeak(W, t.G.null_filtered_white), 1e-5, 'filtered white');
            sd = std(W(:, good), 1, 2);
            t.verifyLessThan(tVortexSignal.relpeak(sd, t.G.null_row_sd), 1e-5, 'row sd');
            % envelope pieces on the stored lag-free output
            sm = @(v) sqrt(rheome.cores.boxsmooth(mean(double(v).^2, 1), 36));
            t.verifyLessThan(tVortexSignal.relpeak(sm(xr), t.G.null_env_data), 1e-5, 'env data');
            t.verifyLessThan(tVortexSignal.relpeak(sm(t.G.null_lagfree), t.G.null_env_lagfree), 1e-5, 'env lagfree');
            ratio = sm(xr) ./ max(sm(t.G.null_lagfree), 1e-30);
            t.verifyLessThan(tVortexSignal.relpeak(ratio, t.G.null_env_ratio), 1e-5, 'env ratio');
        end

        function nullAcceptsRandStream(t)
            xr = single(t.G.frames_real);  good = logical(t.G.good_frames(:).');
            Y = rheome.cores.lagfreenull(xr, good, RandStream('mt19937ar', 'Seed', 1));
            t.verifySize(Y, size(xr));
            t.verifyTrue(all(isfinite(Y), 'all'));
        end
    end
end
