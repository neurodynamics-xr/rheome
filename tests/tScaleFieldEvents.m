classdef tScaleFieldEvents < matlab.unittest.TestCase
% The two MS1 measures end to end on a synthetic cached subject (two icosphere hemispheres, a random
% leadfield, 40 s of noise with a 10 Hz source): fieldsmooth finds the band-limited field smoother than
% the per-vertex one, and eventsensors maps the best Viterbi path onto sensor samples and channels.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_fieldevents'
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            i_cache(fullfile(d, tc.Name));
        end
    end

    methods (Test)

        function bandLimitedFieldIsSmoother(tc)
            S = rheome.scale.sensors(tc.Name);
            [T, X, maps] = rheome.scale.measure_fieldsmooth(tc.Name, S, CutoffsMM=[60 120], NumFrames=4);
            v = @(m, b) T.value(T.metric == m & T.band == b);
            tc.verifyGreaterThan(v("J_wavelength_mm", "bl120"), v("Jt_wavelength_mm", "raw"));
            tc.verifyGreaterThan(v("div_wavelength_mm", "bl120"), v("div_wavelength_mm", "raw"));
            tc.verifyGreaterThan(v("J_coherence", "bl120"), v("Jt_coherence", "raw"));
            tc.verifyGreaterThanOrEqual(v("J_wavelength_mm", "bl120"), v("J_wavelength_mm", "bl60"));
            tc.verifyEqual(height(X), 4 * (7 + 2*5));
            tc.verifySize(maps.bl, [1 2]);  tc.verifySize(maps.raw.Jt, [3*size(maps.surface.Vertices, 1) 1]);
        end

        function eventsLandOnTheirFramesSamples(tc)
            S = rheome.scale.sensors(tc.Name);
            [~, ~, Best] = rheome.scale.measure_grouptrack(tc.Name, S, Depths=3, JLevels=2, Hemis="L");
            tc.assertNotEmpty(Best);
            [T, ev] = rheome.scale.measure_eventsensors(tc.Name, S, Best, Depth=3, JLevel=2);
            tc.verifyEqual(T.value(T.metric == "n_steps"), numel(Best(1).track.frames));
            E = ev.E;  fs = 600;  dec = 2;  bs = 4;
            f1 = (Best(1).window - 1) * 600 + (Best(1).track.frames(1) - 1) * bs;      % 300-fps frames
            s1 = f1 * dec + 1 - round(ev.t(1) * fs);                                      % sample in the crop
            tc.verifyEqual(E.samples(1, :), [s1, s1 + bs*dec - 1]);
            tc.verifyTrue(any(E.mask(:, s1)));  tc.verifyEqual(size(E.mask), size(ev.Xraw));
            tc.verifySize(ev.pos2, [size(ev.Xraw, 1) 2]);
            h = rheome.show.eventsensors(ev.Xband, ev.t, E, ev.pos2, 'Labels', ev.labels, 'Visible', 'off');
            close(h);
        end

        function theDriverWritesTablesAndFigureData(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses=["fieldsmooth" "eventsensors"], OutDir=od);
            tc.verifyEqual(R.timing.status', ["ok" "error"]);           % eventsensors without grouptrack
            tc.verifySubstring(R.timing.message(2), "needs grouptrack first");
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'fieldsmooth.csv')));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'fieldsmooth_maps.mat')));
            tc.verifyTrue(any(R.metrics.metric == "div_wavelength_mm" & R.metrics.band == "bl92"));
        end

        function theDriverWritesPrognomeCoefficients(tc)
            % rheome.scale.run(..., Analyses="coefficients") -> rheome_coeffs.mat, Prognome's contract
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses="coefficients", OutDir=od);
            tc.assertEqual(R.timing.status, "ok", R.timing.message);
            C = load(fullfile(od, tc.Name, 'rheome_coeffs.mat'));
            tc.verifyEmpty(setdiff(["W" "codes" "scales" "fs"], string(fieldnames(C))));
            tc.verifyEqual(C.fs, 100);
            tc.verifyEqual(C.codes, (2^8:2^9-1)');                    % depth 7 per hemisphere: level 8
            tc.verifySize(C.W, [256 numel(C.scales) 40*100]);          % 40 s at 100 Hz
            tc.verifyTrue(all(isfinite(C.W), 'all') && all(C.W >= 0, 'all'));
            tc.verifyEqual(R.metrics.value(R.metrics.metric == "n_tiles"), 256);
        end
    end
end

function i_cache(d)
    mkdir(d);  rng(7);
    [V0, F0] = rheome.geom.icosphere(3);  nh = size(V0, 1);  R = 0.04;
    Vs = {V0*R - [0.03 0 0], V0*R + [0.03 0 0]};
    S = struct('Vertices', [Vs{1}; Vs{2}], 'Faces', [F0; F0 + nh], 'VertNormals', [V0; V0], 'nV', 2*nh, 'nF', 2*size(F0,1));
    bases = struct('hemi', {{'L', 'R'}});
    for h = 1:2
        Sh = struct('Vertices', Vs{h}, 'Faces', F0, 'VertNormals', V0, 'nV', nh, 'nF', size(F0, 1));
        [L, M] = rheome.operators.laplace_beltrami(Sh.Vertices, Sh.Faces);
        [P, D] = eig(full(L), full(M), 'chol');  [lam, o] = sort(max(diag(D), 0));  P = P(:, o(1:200));  lam = lam(1:200);
        P = P ./ sqrt(sum(P .* (M * P), 1));
        lr = 'LR';  bases.(lr(h)) = struct('lbo', struct('L', L, 'M', M, 'Phi', P, 'Lambda', lam, 'Mass', M), ...
                                 'gv', (1:nh)' + (h-1)*nh, 'S', Sh);
    end
    save(fullfile(d, 'bases.mat'), 'bases');  L = [];  M = [];  save(fullfile(d, 'surface.mat'), 'S', 'L', 'M');
    C = 40;  fs = 600;  N = 40 * fs;  az = 2*pi*(0:C-1)/C;  el = linspace(0.2, 1.2, C);
    loc = 0.12 * [cos(el).*cos(az); cos(el).*sin(az); sin(el)];
    Gain = zeros(C, 3*S.nV);                                   % dipole-like falloff, so the gain is smooth
    for c = 1:C
        r = loc(:, c)' - S.Vertices;  g = r ./ vecnorm(r, 2, 2).^3;  Gain(c, :) = reshape(g', 1, []);
    end
    src = 7;  t = (0:N-1)/fs;
    J = zeros(3*S.nV, 1);  J(3*src-2:3*src) = [1 0 0];
    F = 1e-3 * Gain * J * (sin(2*pi*10*t) .* (1 + sin(2*pi*0.2*t))) + 1e-2 * randn(C, N) * std(Gain(:));
    ch = struct('Loc', num2cell(loc, 1)', 'Name', compose('M%02d', (1:C)'));
    chan = struct('Type', {repmat({'MEG'}, C, 1)}, 'Channel', ch);
    rec = struct('F', F, 'sfreq', fs, 'ChannelFlag', ones(C, 1), 'Time', t);
    hm = struct('Gain', Gain);  ncov = struct('NoiseCov', eye(C) * var(F(:)) * 1e-2);
    save(fullfile(d, 'study.mat'), 'rec', 'chan', 'hm', 'ncov');
end

% Author: Diellor Basha, 2026
