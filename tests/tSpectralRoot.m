classdef tSpectralRoot < matlab.unittest.TestCase
% rheome.spectral.root / column / meyerband: the spectral root of a chunk and the column at a time cut from it.
% Synthetic rows only. The method's own checks (the raw column is the sample; the analytic column's imaginary part is the
% Hilbert transform of its real part; neighbouring Meyer bands partition power), then the editor's own functions on the same
% chunk (tests/editor_spectral.mts: fftCpu, bandGains, spectralColumnCpu, spectralAnalyticCpu, compiled from the editor's
% source), skipped when the editor checkout or node are absent.
%
% Author: Diellor Basha, 2026

    properties
        x, halo = 200, rate = 300, band = [8 13], cols = [0 17 1499 2599]
    end

    methods (TestMethodSetup)
        function rows(tc)
            rng = RandStream('mt19937ar', 'Seed', 7); t = (0:2999) / tc.rate;
            tc.x = single(randn(rng, 5, 3000) + 3 * sin(2 * pi * 10 * t + (1:5)') + 40);   % noise, a 10 Hz rhythm, an offset
        end
    end

    methods (Test)
        function rawColumnIsTheSample(tc)
            S = rheome.spectral.root(tc.x, tc.halo);
            tc.verifyEqual([S.N S.width], [4096 3000 - 2 * tc.halo]);
            x0 = rheome.spectral.column(S, tc.cols);
            tc.verifyEqual(x0, double(tc.x(:, tc.halo + tc.cols + 1)), 'AbsTol', 1e-10);
        end

        function analyticImaginaryIsTheHilbertOfTheReal(tc)
            S = rheome.spectral.root(tc.x, tc.halo);
            H = rheome.spectral.meyerband((0:S.N/2)' * tc.rate / S.N, tc.band(1), tc.band(2));
            [xr, xa] = rheome.spectral.column(S, 0:S.N-1-2*tc.halo, H);       % every column of the circular transform's frame
            full = rheome.spectral.column(S, -tc.halo + (0:S.N-1), H);
            h = imag(hilbert(full'))';
            tc.verifyEqual(imag(xa), h(:, tc.halo + 1:S.N - tc.halo), 'AbsTol', 1e-9 * max(abs(h), [], 'all'));
            tc.verifyEqual(real(xa), xr, 'AbsTol', 1e-12);
            amp = median(abs(xa(:, 500:2000)), 2);
            tc.verifyEqual(amp, 3 * ones(5, 1), 'RelTol', 0.15, 'the 10 Hz rhythm''s amplitude, through the band');
        end

        function meyerBandsPartitionPower(tc)
            f = logspace(-1, 2, 400)';
            H1 = rheome.spectral.meyerband(f, 4, 8); H2 = rheome.spectral.meyerband(f, 8, 16);
            on = f > 8 * 2^(-1/3) & f < 8 * 2^(1/3);
            tc.verifyEqual(H1(on).^2 + H2(on).^2, ones(nnz(on), 1), 'AbsTol', 1e-12);
            tc.verifyEqual(rheome.spectral.meyerband([0; 10.2; 1e3], 8, 13), [0; 1; 0], 'AbsTol', 1e-12);   % 10.2 Hz: past 8*2^(1/3) = 10.08, short of 13*2^(-1/3) = 10.32
        end

        function editorFunctionsAgree(tc)
            R = fullfile(getenv('HOME'), 'workspace/research/code/nxr-cortical-flow');
            tc.assumeTrue(isfile(fullfile(R, 'editor/src/editor/engine/gpu/modeView.ts')) && system('command -v npx >/dev/null') == 0, 'no editor checkout or node');
            d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
            fid = fopen(fullfile(d, 'rows.f32'), 'w'); fwrite(fid, tc.x', 'single'); fclose(fid);
            m = struct('C', size(tc.x, 1), 'w', size(tc.x, 2), 'halo', tc.halo, 'rate', tc.rate, 'band', tc.band, 'cols', tc.cols);
            fid = fopen(fullfile(d, 'meta.json'), 'w'); fwrite(fid, jsonencode(m)); fclose(fid);
            [st, out] = system(sprintf('cd "%s/editor" && npx --no-install tsx "%s" "%s" "%s" 2>&1', R, fullfile(fileparts(mfilename('fullpath')), 'editor_spectral.mts'), R, d));
            tc.assertEqual(st, 0, out);
            S = rheome.spectral.root(tc.x, tc.halo); C = m.C; bins = S.N / 2 + 1;
            Xe = reshape(i_read(fullfile(d, 'X.f32'), 'single'), 2, bins, C); Xe = squeeze(complex(Xe(1,:,:), Xe(2,:,:))).';
            He = i_read(fullfile(d, 'H.f64'), 'double');
            E = reshape(i_read(fullfile(d, 'cols.f64'), 'double'), 3, C, []);
            tc.verifyLessThan(max(abs(S.X - Xe), [], 'all') / max(abs(Xe), [], 'all'), 1e-7, 'the root: the editor holds it in single');
            H = rheome.spectral.meyerband((0:S.N/2)' * tc.rate / S.N, tc.band(1), tc.band(2));
            tc.verifyEqual(H, He, 'AbsTol', 1e-7, 'bandGains (single)');
            S.X = Xe; x0 = rheome.spectral.column(S, tc.cols); [~, xa] = rheome.spectral.column(S, tc.cols, He);   % the editor's single inputs
            sc = @(a, b) max(abs(a - b), [], 'all') / max(abs(b), [], 'all');
            tc.verifyLessThan(sc(x0, squeeze(E(1,:,:))), 1e-6, 'spectralColumnCpu');
            tc.verifyLessThan(sc(real(xa), squeeze(E(2,:,:))), 1e-6, 'spectralAnalyticCpu, real');
            tc.verifyLessThan(sc(imag(xa), squeeze(E(3,:,:))), 1e-6, 'spectralAnalyticCpu, imaginary');
        end
    end
end

function x = i_read(p, t)
fid = fopen(p); x = fread(fid, inf, [t '=>double']); fclose(fid);
end
