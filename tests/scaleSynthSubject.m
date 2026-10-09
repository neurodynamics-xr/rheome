function scaleSynthSubject(d)
% SCALESYNTHSUBJECT  A synthetic cached subject for the rheome.scale measures: two icosphere hemispheres
% (radius 40 mm), a dipole-falloff leadfield on 40 channels, 40 s of noise with a 10 Hz source, and an
% empty room (noise.mat) on the same channels, and a four-region 'Desikan-Killiany' atlas (atlas.mat: each
% hemisphere split at its mean z).
%
%   scaleSynthSubject(fullfile(getenv('RHEOME_DATA'), name))
%
% Shared by tScaleFieldEvents, tScalePlantFloors and tScaleGroupP2.
%
% Author: Diellor Basha, 2026

    mkdir(d);  rng(7);
    [V0, F0] = rheome.geom.icosphere(3);  nh = size(V0, 1);  R = 0.04;
    Vs = {V0*R - [0.03 0 0], V0*R + [0.03 0 0]};
    S = struct('Vertices', [Vs{1}; Vs{2}], 'Faces', [F0; F0 + nh], 'VertNormals', [V0; V0], 'nV', 2*nh, 'nF', 2*size(F0,1), ...
               'Hemi', {{(1:nh)', (nh+1:2*nh)'}}, 'HemiLabel', {{'Cortex L', 'Cortex R'}}, 'Comment', 'synth', 'SurfaceFile', '', ...
               'Sphere', [V0; V0]*0.1);                   % the registration sphere: each hemisphere on R = 100 mm
    bases = struct('hemi', {{'L', 'R'}});
    for h = 1:2
        Sh = struct('Vertices', Vs{h}, 'Faces', F0, 'VertNormals', V0, 'nV', nh, 'nF', size(F0, 1));
        [L, M] = rheome.operators.laplace_beltrami(Sh.Vertices, Sh.Faces);
        [P, D] = eig(full(L), full(M), 'chol');  [lam, o] = sort(max(diag(D), 0));  P = P(:, o(1:200));  lam = lam(1:200);
        P = P ./ sqrt(sum(P .* (M * P), 1));
        lr = 'LR';  bases.(lr(h)) = struct('lbo', struct('L', L, 'M', M, 'Phi', P, 'Lambda', lam, 'Mass', M), ...
                                 'gv', (1:nh)' + (h-1)*nh, 'S', Sh);
    end
    save(fullfile(d, 'bases.mat'), 'bases');  [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces);  save(fullfile(d, 'surface.mat'), 'S', 'L', 'M');
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
    chan = struct('Type', {repmat({'MEG'}, 1, C)}, 'Name', {{ch.Name}}, 'Channel', ch');
    rec = struct('F', F, 'sfreq', fs, 'ChannelFlag', ones(C, 1), 'Time', t);
    hm = struct('Gain', Gain, 'nV', S.nV);  ncov = struct('NoiseCov', eye(C) * var(F(:)) * 1e-2, 'FourthMoment', [], 'nSamples', []);
    save(fullfile(d, 'study.mat'), 'rec', 'chan', 'hm', 'ncov');
    nrec = struct('F', 1e-2 * randn(C, N) * std(Gain(:)), 'sfreq', fs, 'ChannelFlag', ones(C, 1), 'Time', t, ...
                  'ChannelName', {{ch.Name}'}, 'ChannelType', {chan.Type});
    save(fullfile(d, 'noise.mat'), 'nrec');
    z = S.Vertices(:, 3) > mean(S.Vertices(:, 3));  hL = (1:S.nV)' <= nh;
    atlas = struct('Name', 'Desikan-Killiany', 'Label', {{'a L', 'b L', 'a R', 'b R'}}, ...
                   'Membership', sparse([z & hL, ~z & hL, z & ~hL, ~z & ~hL]'));
    save(fullfile(d, 'atlas.mat'), 'atlas');
end

% Author: Diellor Basha, 2026
