classdef tVortexTiles < matlab.unittest.TestCase
    % rheome.cores surface ladder, tile crossings, path consistency R, repeat, tile dot, tile p/q, BH and the place-based direction
    % consistency against nxr's synthetic goldens (tests/goldens/vortex-tiles, nxr-analysis commit 0db82e4). nxr is the reference.
    %
    % Permutations are replayed from the stored ones (stage_perms, dc_perms): nxr's numpy generator cannot be reproduced.
    %
    % Author: Diellor Basha, 2026

    properties
        G
        Labels   % [V x 9]
        Perms    % cell(3 levels, 4 stages) of cell(5,1) 1-based permutations
        Levels = [6 7 8]
        Stages = {'burst', 'onset', 'within', 'offset'}
    end

    methods (TestClassSetup)
        function load(t)
            d = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'vortex-tiles');
            M = jsondecode(fileread(fullfile(d, 'manifest.json')));
            nm = fieldnames(M.arrays);
            for i = 1:numel(nm), t.G.(nm{i}) = tVortexTiles.readArray(d, M.arrays.(nm{i})); end
            t.Labels = double(t.G.labels).';
            len = double(t.G.stage_perm_lengths(:));  p = double(t.G.stage_perms(:));  off = [0; cumsum(len)];
            t.Perms = cell(3, 4);  c = 0;
            for l = 1:3
                for s = 1:4
                    t.Perms{l, s} = cell(5, 1);
                    for k = 1:5, c = c + 1;  t.Perms{l, s}{k} = p(off(c) + 1 : off(c + 1)).' + 1; end
                end
            end
        end
    end

    methods (Static)
        function x = readArray(d, a)
            sh = double(a.shape(:)).';
            f = fopen(fullfile(d, a.file), 'r'); x = fread(f, prod(sh), ['*' strrep(a.dtype, 'float64', 'double')]); fclose(f);
            if numel(sh) == 1, x = x(:); else, x = permute(reshape(x, fliplr(sh)), numel(sh):-1:1); end
        end
    end

    methods
        function [sums, n, key, vec] = cross(t, S, frames, l)
            [sums, n, key, vec] = rheome.cores.crossings(S, frames, t.Labels(:, l + 1), t.G.V_mm, t.G.normals, double(t.G.support) + 1);
        end
    end

    % The stored per-run sums are float32 (as nxr keeps them) and nxr's repeat / tile dot / R on them run in float32; here they run in double,
    % so those outputs agree to float32 rounding (1e-6), not 1e-12. Crossings, R, chance and place consistency, computed from the steps, agree to 1e-12.
    methods (Test)
        function ladderCodesExact(t)
            lab = rheome.cores.ladder(t.G.V_mm, double(t.G.F) + 1, t.G.hemi);
            t.verifyClass(lab, 'int32');
            t.verifySize(lab, [1284 9]);
            t.verifyEqual(double(lab), t.Labels, 'ladder codes');
            for l = 1:8, t.verifyEqual(floor(double(lab(:, l + 1)) / 2), double(lab(:, l))); end
            t.verifyEqual(numel(unique(lab(:, 9))), double(t.G.tile_count_per_level(9)));
        end

        function ladderHemisphereLetters(t)
            h = repmat("L", 1284, 1);  h(t.G.hemi == 1) = "R";
            lab = rheome.cores.ladder(t.G.V_mm, double(t.G.F) + 1, h);
            t.verifyEqual(double(lab), t.Labels);
        end

        function nearestSupportVertex(t)
            S = double(single(t.G.steps_A));  Vs = t.G.V_mm(double(t.G.support) + 1, :);
            [~, i] = min(sum((permute(S(:, 2:4), [1 3 2]) - permute(Vs, [3 1 2])).^2, 3), [], 2);
            t.verifyEqual(i - 1, double(t.G.nearest_support_A));
        end

        function crossingsExact(t)
            fr = double(t.G.frames_burst);
            for l = t.Levels
                [s, n, k, v] = t.cross(t.G.steps_A, fr, l);
                t.verifyEqual(n, double(t.G.(sprintf('cross_n_A_l%d', l))), sprintf('n l%d', l));
                t.verifyLessThan(max(abs(s - t.G.(sprintf('cross_sums_A_l%d', l))), [], 'all'), 1e-12, sprintf('sums l%d', l));
                t.verifyEqual(n, double(t.G.(sprintf('n_A_l%d', l))));
                if l ~= 7
                    t.verifyEqual(k - 1, double(t.G.(sprintf('cross_key_A_l%d', l))), sprintf('order l%d', l));
                    t.verifyLessThan(max(abs(v - t.G.(sprintf('cross_vec_A_l%d', l))), [], 'all'), 1e-12, sprintf('vec l%d', l));
                end
            end
        end

        function pathRAndChance(t)
            R = zeros(3, 4);  sh = zeros(3, 4);  nc = zeros(3, 4);
            for li = 1:3
                l = t.Levels(li);
                for si = 1:4
                    [s, n, k, v] = t.cross(t.G.steps_A, double(t.G.(['frames_' t.Stages{si}])), l);
                    R(li, si) = rheome.cores.pathR(s, n);  nc(li, si) = sum(n);
                    sh(li, si) = rheome.cores.pathchance(k, v, numel(n), t.Perms{li, si});
                end
            end
            t.verifyEqual(nc, double(t.G.stage_crossings_A), 'crossings');
            t.verifyLessThan(max(abs(R - t.G.stage_R_A), [], 'all'), 1e-12, 'R');
            t.verifyLessThan(max(abs(sh - t.G.stage_shuf_A), [], 'all'), 1e-12, 'chance');
            disp(table(max(abs(R - t.G.stage_R_A), [], 'all'), max(abs(sh - t.G.stage_shuf_A), [], 'all'), 'VariableNames', {'maxErrR', 'maxErrChance'}));
        end

        function burstRofThreeRuns(t)
            for r = 'ABC'
                for li = 1:3
                    l = t.Levels(li);
                    R = rheome.cores.pathR(t.G.(sprintf('sums_%s_l%d', r, l)), double(t.G.(sprintf('n_%s_l%d', r, l))));
                    % nxr stores the float32 sums; its R was taken from the float64 ones (<= 7e-10 apart)
                    t.verifyLessThan(abs(R - t.G.(['burst_R_' r])(li)), 1e-8, sprintf('burst R %s l%d', r, l));
                end
            end
            [s, n] = t.crossB('B', 6);
            t.verifyLessThan(max(abs(s - t.G.sums_B_l6), [], 'all'), 5e-6);   % the stored sums are float32
            t.verifyEqual(n, double(t.G.n_B_l6));
        end

        function repeatTileDotPQ(t)
            for l = t.Levels
                A = {t.G.(sprintf('sums_A_l%d', l)), double(t.G.(sprintf('n_A_l%d', l)))};
                B = {t.G.(sprintf('sums_B_l%d', l)), double(t.G.(sprintf('n_B_l%d', l)))};
                C = {t.G.(sprintf('sums_C_l%d', l)), double(t.G.(sprintf('n_C_l%d', l)))};
                t.verifyLessThan(abs(rheome.cores.repeat(A{:}, B{:}) - t.G.(sprintf('repeat_AB_l%d', l))), 1e-6, 'repeat AB');
                t.verifyLessThan(abs(rheome.cores.repeat(A{:}, C{:}) - t.G.(sprintf('repeat_AC_l%d', l))), 1e-6, 'repeat AC');
                dab = rheome.cores.tiledot(A{:}, B{:});  dac = rheome.cores.tiledot(A{:}, C{:});
                t.verifyEqual(isnan(dab), isnan(t.G.(sprintf('tile_dot_AB_l%d', l))), 'NaN pattern AB');
                t.verifyLessThan(max(abs(dab - t.G.(sprintf('tile_dot_AB_l%d', l))), [], 'omitnan'), 1e-6, 'dot AB');
                t.verifyEqual(isnan(dac), isnan(t.G.(sprintf('tile_dot_AC_l%d', l))));
                t.verifyLessThan(max(abs(dac - t.G.(sprintf('tile_dot_AC_l%d', l))), [], 'omitnan'), 1e-6, 'dot AC');
                [p, q] = rheome.cores.tilepq(dab, t.G.(sprintf('tile_dot_null_l%d', l)));
                pr = t.G.(sprintf('tile_p_AB_l%d', l));  qr = t.G.(sprintf('tile_q_AB_l%d', l));
                t.verifyEqual(isnan(p), isnan(pr));  t.verifyEqual(isnan(q), isnan(qr));
                t.verifyLessThan(max(abs(p - pr), [], 'omitnan'), 1e-12, 'p');
                t.verifyLessThan(max(abs(q - qr), [], 'omitnan'), 1e-12, 'q');
            end
        end

        function benjaminiHochberg(t)
            q = rheome.cores.bh(t.G.bh_p);
            t.verifyLessThan(max(abs(q - t.G.bh_q)), 1e-12);
            t.verifyEqual(q(:).', [0.05 0.0667 0.0667 0.2857 0.625 0.0667 1 0.01 0.0667 1], AbsTol=5e-5);
        end

        function placeConsistency(t)
            S = t.G.steps_A(t.G.steps_A(:, 8) < 2000, :);
            dp = double(t.G.dc_perms(:)) ;  len = double(t.G.dc_perm_lengths(:));  off = [0; cumsum(len)];
            perms = arrayfun(@(i) dp(off(i) + 1 : off(i + 1)).' + 1, 1:numel(len), 'UniformOutput', false);
            [R, ch, nu] = rheome.cores.placeconsistency(S, t.G.V_mm, t.G.normals, double(t.G.dc_places) + 1, perms, 10);
            t.verifyEqual(nu, t.G.dc_result(3));
            t.verifyLessThan(abs(R - t.G.dc_result(1)), 1e-12, 'R');
            t.verifyLessThan(abs(ch - t.G.dc_result(2)), 1e-12, 'chance');
            t.verifyGreaterThan(R, ch + 0.2);
        end
    end

    methods
        function [s, n] = crossB(t, run, l)
            [s, n] = t.cross(t.G.(['steps_' run]), double(t.G.frames_burst), l);
        end
    end
end
