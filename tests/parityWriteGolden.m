function out = parityWriteGolden(fixture, out)
% PARITYWRITEGOLDEN  Rheome's numbers on the Cortical Flow Editor's cortex fixture, for the Editor's parity test (rheome.golden@1).
%
%   parityWriteGolden(fixture, out)
%     fixture : the Editor's fixture stem, e.g. '<nxr-cortical-flow>/editor/fixtures/cortex_pial_low'
%               (reads <stem>.positions.f32, float32 [3 x nV], and <stem>.faces.u32, uint32 0-based)
%     out     : the golden folder (created). The Editor reads it: RHEOME_GOLDEN=<out> in editor/src/parity/rheome.parity.test.ts
%
% A Gaussian bump (sigma 12 mm) moves 5 mm a frame along +y for 21 frames at 100 Hz. Written as float64 little-endian,
% ROW-MAJOR [rows x cols] as golden.json names them: Lx (mesh cotan L * frame 1), Mlumped, MgalRowSum, I [T x nV],
% flow_alpha1 and flow_alpha0_01 [(T-1) x 3nV] (rheome.dynamics.opticalflow_scalar), div and curl of flow_alpha1
% [(T-1) x nV] (the face gradient). Not a test: the parity test lives on the Editor's side; Rheome is the reference.
%
% Author: Diellor Basha, 2026

    if ~isfolder(out), mkdir(out); end
    fid = fopen([fixture '.positions.f32'], 'r', 'ieee-le'); V = double(reshape(fread(fid, inf, 'single'), 3, [])'); fclose(fid);
    fid = fopen([fixture '.faces.u32'], 'r', 'ieee-le'); F = double(reshape(fread(fid, inf, 'uint32'), 3, [])') + 1; fclose(fid);
    fs = 100; T = 21; c0 = [0.04, -0.05, 0.06]; step = [0, 0.005, 0]; sigma = 0.012;
    I = zeros(size(V, 1), T);
    for t = 1:T, I(:, t) = exp(-sum((V - (c0 + (t-1)*step)).^2, 2) / (2*sigma^2)); end
    s = struct('Vertices', V, 'Faces', F, 'VertNormals', vertexNormal(triangulation(F, V)), 'nV', size(V, 1));

    [L, Mg] = rheome.operators.laplace_beltrami(V, F, 'galerkin');
    [~, Ml] = rheome.operators.laplace_beltrami(V, F, 'lumped');
    G = struct();
    G.Lx = put(out, 'Lx', L*I(:, 1));
    G.Mlumped = put(out, 'Mlumped', full(diag(Ml)));
    G.MgalRowSum = put(out, 'MgalRowSum', full(sum(Mg, 2)));
    G.I = put(out, 'I', I');
    tic; v1 = rheome.dynamics.opticalflow_scalar(I, s, struct('alpha', 1));    t1 = toc;
    tic; v2 = rheome.dynamics.opticalflow_scalar(I, s, struct('alpha', 0.01)); t2 = toc;
    G.flow_alpha1 = put(out, 'flow_alpha1', v1');
    G.flow_alpha0_01 = put(out, 'flow_alpha0_01', v2');
    fg = rheome.operators.face_gradient(V, F);
    G.div = put(out, 'div', rheome.differential.divergence(v1, s, fg)');
    G.curl = put(out, 'curl', rheome.differential.curl(v1, s, fg)');

    r = regexp(fileread(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Contents.m')), '% Version (\S+)', 'tokens', 'once');   % not ver(): the clone's folder name may differ
    m = struct('format', 'rheome.golden@1', 'fixture', 'editor/fixtures/cortex_pial_low', 'nV', size(V, 1), 'T', T, 'fs', fs, ...
        'bump', struct('c0', c0, 'step_per_frame', step, 'sigma', sigma), 'rheome', r{1}, 'matlab', version('-release'), ...
        'seconds', struct('flow_alpha1', t1, 'flow_alpha0_01', t2), 'arrays', G);
    fid = fopen(fullfile(out, 'golden.json'), 'w'); fwrite(fid, jsonencode(m, 'PrettyPrint', true)); fclose(fid);
end

function d = put(out, name, A)
    fid = fopen(fullfile(out, [name '.bin']), 'w', 'ieee-le'); fwrite(fid, A', 'double'); fclose(fid);   % row-major
    d = struct('file', [name '.bin'], 'shape', size(A));
end
