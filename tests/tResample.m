classdef tResample < matlab.unittest.TestCase
% rheome.filters.resample: the length rule, a band-limited sine carried through, the class kept, and scipy's resample_poly
% (tests/resample_poly.py, nxr-convert's call) on synthetic rows for several ratios, down and up. The scipy check is skipped
% when nxr-convert or uv are absent.
%
% Author: Diellor Basha, 2026

    properties (TestParameter)
        rates = {[1200 300], [300 1200], [1000 300], [2400 600], [22050 300], [250 200]}
    end

    methods (Test)
        function lengthAndRatio(tc)
            [y, p, q] = rheome.filters.resample(zeros(2, 1001), 1200, 300);
            tc.verifyEqual([size(y) p q], [2 251 1 4]);
            tc.verifyError(@() rheome.filters.resample(zeros(1, 10), 1000, 999.9), 'rheome:resample');
            tc.verifyEqual(rheome.filters.resample(ones(3, 5), 300, 300), ones(3, 5));
            tc.verifyClass(rheome.filters.resample(single(ones(2, 100)), 1200, 300), 'single');
            tc.verifySize(rheome.filters.resample(ones(100, 2), 1200, 300, Dim=1), [25 2]);
        end

        function sineIsCarried(tc)
            t = (0:11999) / 1200; x = sin(2 * pi * 10 * t + [0; 1]);
            y = rheome.filters.resample(x, 1200, 300);
            u = (0:size(y, 2) - 1) / 300; ref = sin(2 * pi * 10 * u + [0; 1]);
            mid = 50:size(y, 2) - 50;                                          % away from the zero-padded ends
            tc.verifyLessThan(max(abs(y(:, mid) - ref(:, mid)), [], 'all'), 1e-3);
        end

        function scipyAgrees(tc, rates)
            nxr = fullfile(getenv('HOME'), 'workspace/research/code/nxr-cortical-flow/backend/convert/nxr-convert');
            tc.assumeTrue(isfile(fullfile(nxr, 'pyproject.toml')) && system('command -v uv >/dev/null') == 0, 'nxr-convert or uv not found');
            x = randn(RandStream('mt19937ar', 'Seed', 11), 3, 1001) + 5;
            d = tempname; mkdir(d); c = onCleanup(@() rmdir(d, 's'));
            fid = fopen(fullfile(d, 'x.f64'), 'w'); fwrite(fid, x', 'double'); fclose(fid);
            [st, out] = system(sprintf('uv run --frozen --project "%s" python -I "%s" "%s" 3 %g %g "%s" 2>&1', nxr, ...
                fullfile(fileparts(mfilename('fullpath')), 'resample_poly.py'), fullfile(d, 'x.f64'), rates(1), rates(2), fullfile(d, 'y.f64')));
            tc.assertEqual(st, 0, out);
            fid = fopen(fullfile(d, 'y.f64')); ys = fread(fid, inf, 'double'); fclose(fid);
            y = rheome.filters.resample(x, rates(1), rates(2));
            tc.assertEqual(numel(y), numel(ys), 'the same length');
            tc.verifyEqual(y, reshape(ys, [], 3)', 'AbsTol', 1e-12 * max(abs(ys)));
        end
    end
end
