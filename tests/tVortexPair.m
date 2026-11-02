classdef tVortexPair < matlab.unittest.TestCase
    % rheome.cores lagged pair, share, S3 and Wishart null against nxr's synthetic goldens (tests/goldens/vortex-pair,
    % nxr-analysis commit 3da1320; step 10a). nxr is the reference. The eigenvector phase is free, so raw fa, fb, u are
    % compared in the canonical gauge (rheome.cores.pairgauge); gauge-free invariants are compared directly. Stored Wishart
    % draws are replayed.
    %
    % Author: Diellor Basha, 2026

    properties
        G
        P
        Scn = {'rot', 'cw', 'ctl'}
        Op = {'I', 'L'}
    end

    methods (TestClassSetup)
        function load(t)
            d = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'vortex-pair');
            M = jsondecode(fileread(fullfile(d, 'manifest.json')));
            nm = fieldnames(M.arrays);
            for i = 1:numel(nm), t.G.(nm{i}) = tVortexPair.readArray(d, M.arrays.(nm{i})); end
            t.P = M.parameters;
        end
    end

    methods (Static)
        function x = readArray(d, a)
            sh = double(a.shape(:)).';
            f = fopen(fullfile(d, a.file), 'r'); x = fread(f, max(prod(sh), 1), ['*' strrep(a.dtype, 'float64', 'double')]); fclose(f);
            if numel(sh) == 0, x = double(x); elseif numel(sh) == 1, x = x(:); else, x = permute(reshape(x, fliplr(sh)), numel(sh):-1:1); end
        end

        function e = relpeak(x, ref)
            e = max(abs(x(:) - ref(:))) / max(abs(ref(:)));
        end
    end

    methods
        function S = cx(t, s, name), S = t.G.([s '_' name '_real']) + 1i * t.G.([s '_' name '_imag']); end
        function L = op(t, o), L = t.G.(['operator_' o]); end
        function [fa, fb, u] = fit(t, s, o)
            [fa, fb, u] = rheome.cores.pair(t.op(o), t.cx(s, 'S_even'));
        end
    end

    methods (Test)
        function covariances(t)
            for s = t.Scn
                xa = t.cx(s{1}, 'analytic');
                for nm = {'all', 'even', 'odd'}
                    m = {t.G.([s{1} '_ok']), t.G.([s{1} '_even']), t.G.([s{1} '_odd'])};
                    mk = m{strcmp(nm{1}, {'all', 'even', 'odd'})};
                    S = rheome.cores.paircov(xa, mk);
                    t.verifyLessThan(tVortexPair.relpeak(S, t.cx(s{1}, ['S_' nm{1}])), 1e-12, [s{1} ' S_' nm{1}]);
                end
            end
        end

        function dof(t)
            for s = t.Scn
                xa = t.cx(s{1}, 'analytic'); ok = logical(t.G.([s{1} '_ok']));
                f = rheome.cores.pairdof(xa(:, ok));
                t.verifyEqual(f, t.G.([s{1} '_dof_fraction']), 'RelTol', 1e-10);
                od = logical(t.G.([s{1} '_odd']));
                t.verifyEqual([nnz(ok) nnz(od)] * f, t.G.([s{1} '_dof']).', 'RelTol', 1e-10);
            end
        end

        function sqrtPsd(t)
            for s = t.Scn
                R = rheome.cores.pairsqrt(real(t.cx(s{1}, 'S_all')));
                t.verifyLessThan(tVortexPair.relpeak(R, t.G.([s{1} '_sqrt_ReS'])), 1e-10, s{1});
            end
        end

        function analyticSignal(t)
            % no stored sensor input: a planted band-limited sinusoid gives A exp(i phase); equals bandframes at 100 fps
            fs = t.P.fs_hz; n = 3000; tt = (0:n-1) / fs;
            X = [1.5 * cos(2*pi*10*tt + 0.3); 0.7 * sin(2*pi*10*tt)];
            z = rheome.cores.analytic(X, fs, t.P.f0_hz, t.P.fwhm_hz, 3);
            t.verifySize(z, [2 1000]);
            ref = [1.5 * exp(1i*(2*pi*10*tt(1:3:end) + 0.3)); 0.7 * exp(1i*(2*pi*10*tt(1:3:end) - pi/2))];
            t.verifyLessThan(max(abs(z(:) - ref(:))), 1e-10);
            zf = rheome.cores.bandframes(X, fs, t.P.f0_hz, t.P.fwhm_hz, 100);
            t.verifyLessThan(max(abs(z(:) - double(zf(:)))), 1e-5);
            % the stored analytic signal has the planted band's covariance structure
            t.verifyClass(z, 'double');
        end

        function pairRawInGauge(t)
            for s = t.Scn(1:2)
                for o = t.Op
                    [fa, fb, u] = t.fit(s{1}, o{1});
                    [fa, fb, u] = rheome.cores.pairgauge(fa, fb, u);
                    pre = [s{1} '_' o{1} '_'];
                    pk = max(abs(t.G.([pre 'fa'])));
                    t.verifyLessThan(max(abs(fa(:) - t.G.([pre 'fa']))) / pk, 1e-6, [pre 'fa']);
                    t.verifyLessThan(max(abs(fb(:) - t.G.([pre 'fb']))) / pk, 1e-6, [pre 'fb']);
                    uref = t.G.([pre 'u_real']) + 1i * t.G.([pre 'u_imag']);
                    t.verifyLessThan(max(abs(u(:) - uref)) / max(abs(uref)), 1e-6, [pre 'u']);
                end
            end
        end

        function pairInvariants(t)
            for s = t.Scn
                for o = t.Op
                    [fa, fb] = t.fit(s{1}, o{1}); pre = [s{1} '_' o{1} '_'];
                    proj = fa.' * fa + fb.' * fb; form = fa.' * fb - fb.' * fa;
                    sc = max(abs(t.G.([pre 'proj'])), [], 'all');
                    if ~strcmp(s{1}, 'ctl')       % in ctl the pair is a noise artefact (small gap)
                        t.verifyLessThan(max(abs(proj - t.G.([pre 'proj']))/sc, [], 'all'), 1e-6, [pre 'proj']);
                        t.verifyLessThan(max(abs(form - t.G.([pre 'form']))/sc, [], 'all'), 1e-6, [pre 'form']);
                    end
                    t.verifyGreaterThanOrEqual(imag(fa * t.cx(s{1}, 'S_even') * fb.'), 0, [pre 'orientation']);
                end
            end
        end

        function shareAndS3(t)
            for s = t.Scn
                for o = t.Op
                    pre = [s{1} '_' o{1} '_']; L = t.op(o{1});
                    t.verifyEqual(rheome.cores.pairshare(L, t.cx(s{1}, 'S_all')), t.G.([pre 'share']), 'AbsTol', 1e-10, [pre 'share']);
                    t.verifyEqual(rheome.cores.pairshare(L, t.cx(s{1}, 'S_even')), t.G.([pre 'share_even']), 'AbsTol', 1e-10, [pre 'share_even']);
                    t.verifyEqual(rheome.cores.pairshare(L, t.cx(s{1}, 'S_pop')), t.G.([pre 'share_pop']), 'AbsTol', 1e-10, [pre 'share_pop']);
                    [fa, fb] = t.fit(s{1}, o{1});
                    So = t.cx(s{1}, 'S_odd');
                    t.verifyEqual(rheome.cores.pairs3(fa, fb, So), t.G.([pre 's3']), 'AbsTol', 1e-8, [pre 's3']);
                    t.verifyEqual(rheome.cores.pairs3(fa, fb, t.cx(s{1}, 'S_even')), t.G.([pre 's3_even']), 'AbsTol', 1e-8, [pre 's3_even']);
                end
            end
        end

        function wishartReplay(t)
            for s = t.Scn
                Rh = t.G.([s{1} '_sqrt_ReS']); dof = t.G.([s{1} '_dof']);
                Da = t.G.([s{1} '_draws_all']); Do = t.G.([s{1} '_draws_odd']);   % [N 2 C k], as stored
                N = size(Da, 1);
                t.assertEqual(N, t.P.nnull);
                NA = t.cx(s{1}, 'null_S_all'); NO = t.cx(s{1}, 'null_S_odd');
                for i = [1 2 17 N]
                    da = permute(squeeze(Da(i, :, :, :)), [2 3 1]); do = permute(squeeze(Do(i, :, :, :)), [2 3 1]);
                    Sa = rheome.cores.pairnull(Rh, dof(1), da); So = rheome.cores.pairnull(Rh, dof(2), do);
                    t.verifyLessThan(tVortexPair.relpeak(Sa, squeeze(NA(i, :, :))), 1e-10, 'null all');
                    t.verifyLessThan(tVortexPair.relpeak(So, squeeze(NO(i, :, :))), 1e-10, 'null odd');
                end
            end
        end

        function pValues(t)
            for s = t.Scn
                Rh = t.G.([s{1} '_sqrt_ReS']); dof = t.G.([s{1} '_dof']);
                Da = t.G.([s{1} '_draws_all']); Do = t.G.([s{1} '_draws_odd']); N = size(Da, 1);
                for o = t.Op
                    pre = [s{1} '_' o{1} '_']; L = t.op(o{1});
                    [fa, fb] = t.fit(s{1}, o{1});
                    shn = zeros(N, 1); s3n = zeros(N, 1);
                    for i = 1:N
                        da = permute(squeeze(Da(i, :, :, :)), [2 3 1]); do = permute(squeeze(Do(i, :, :, :)), [2 3 1]);
                        shn(i) = rheome.cores.pairshare(L, rheome.cores.pairnull(Rh, dof(1), da));
                        s3n(i) = rheome.cores.pairs3(fa, fb, rheome.cores.pairnull(Rh, dof(2), do));
                    end
                    t.verifyEqual(shn, t.G.([pre 'share_null']), 'AbsTol', 1e-9, [pre 'share null']);
                    t.verifyEqual(s3n, t.G.([pre 's3_null']), 'AbsTol', 1e-8, [pre 's3 null']);
                    p1 = (1 + nnz(shn >= t.G.([pre 'share']))) / (N + 1);
                    p2 = (1 + nnz(s3n >= t.G.([pre 's3']))) / (N + 1);
                    t.verifyEqual(p1, t.G.([pre 'p_share']), 'AbsTol', 1e-12, [pre 'p_share']);
                    t.verifyEqual(p2, t.G.([pre 'p_s3']), 'AbsTol', 1e-12, [pre 'p_s3']);
                end
            end
        end

        function plantedAnswers(t)
            p = t.G.planted_p(:); q = t.G.planted_q(:);
            for s = t.Scn(1:2)
                [fa, fb] = t.fit(s{1}, 'I');
                A = [fa(:) fb(:)]; B = [p q];
                Qa = orth(A); Qb = orth(B);
                ang = acosd(min(1, min(svd(Qa' * Qb))));
                t.verifyLessThan(ang, 8, [s{1} ' plane angle']);
                t.verifyGreaterThan(t.G.([s{1} '_I_s3']), 0.5, [s{1} ' held-out S3']);
                t.verifyEqual(t.G.([s{1} '_I_p_share']), 1 / (t.P.nnull + 1), 'AbsTol', 1e-12);
            end
            % the control (standing pattern): no lag, small share, nothing significant
            t.verifyLessThan(t.G.ctl_I_share, 0.2);
            t.verifyGreaterThan(t.G.ctl_I_p_share, 0.05);
            % a planted quarter-cycle pair is recovered exactly from its noise-free cross-spectrum
            C = 5; a = [1 0 0 0 0]; b = [0 1 0 0 0];
            [fa, fb] = rheome.cores.pair(eye(C), a.' * a + b.' * b + 1i * (b.' * a - a.' * b));
            t.verifyEqual(fa.' * fa + fb.' * fb, diag([0.5 0.5 0 0 0]), 'AbsTol', 1e-12);
            t.verifyEqual(rheome.cores.pairs3(fa, fb, a.' * a + b.' * b + 1i * (b.' * a - a.' * b)), 1, 'AbsTol', 1e-12);
        end
    end
end
