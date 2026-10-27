classdef tNxrStoreParity < matlab.unittest.TestCase
% rheome.io.read.nxrsurface -> rheome.geom.glue -> rheome.operators.hodge_dirac -> rheome.eigen.hodge_dirac against nxr-convert's
% own chain (tests/nxr_convert_chain.py: modes_input, glue_mesh, prepared, solve_hodge) on every subject of a real nxr store:
% the same surface and wall read, the same glued mesh (bit for bit) and key, the same flips, the same band of modes.
%
% The store is named by $RHEOME_NXR_STORE (a dataset folder holding dataset.sqlite) and never written: nxr-convert reads it through
% a read-only sqlite URI. Nothing from it is kept: no subject, value or array of a store lives in this repository (it is public).
% Skipped when the variable is unset, the store is absent, or nxr-convert / uv are not found.
%
% Author: Diellor Basha, 2026

    properties
        store, nxr, subjects
    end

    methods (TestClassSetup)
        function find(tc)
            tc.store = string(getenv('RHEOME_NXR_STORE'));
            tc.assumeTrue(tc.store ~= "" && isfile(fullfile(tc.store, 'dataset.sqlite')), 'RHEOME_NXR_STORE is unset or holds no nxr store');
            tc.nxr = fullfile(getenv('HOME'), 'workspace/research/code/nxr-cortical-flow/backend/convert/nxr-convert');
            tc.assumeTrue(isfile(fullfile(tc.nxr, 'pyproject.toml')) && system('command -v uv >/dev/null') == 0, 'nxr-convert or uv not found');
            [st, out] = system(sprintf('sqlite3 -readonly "%s" "SELECT name FROM subject WHERE kind = ''subject'' AND status = ''complete'' ORDER BY name"', ...
                fullfile(tc.store, 'dataset.sqlite')));
            tc.assertEqual(st, 0, out);
            tc.subjects = splitlines(strtrim(string(out)));
            tc.assumeNotEmpty(tc.subjects(tc.subjects ~= ""), 'the store has no complete subject');
        end
    end

    methods (Test)
        function chainIsNxrConverts(tc)
            for k = 1:numel(tc.subjects)
                s = tc.subjects(k); at = sprintf('subject %d', k);           % by number: the diagnostics name no participant
                o = tempname; c = onCleanup(@() rmdir(o, 's'));
                cmd = sprintf('uv run --frozen --project "%s" python -I "%s" "%s" "%s" 200 "%s" 2>&1', tc.nxr, ...
                    fullfile(fileparts(mfilename('fullpath')), 'nxr_convert_chain.py'), tc.store, s, o);
                [st, out] = system(cmd);
                tc.assertEqual(st, 0, sprintf('%s: nxr-convert failed: %s', at, out));
                py = @(n, w) reshape(i_read(fullfile(o, n + ".f64")), w, [])';

                S = rheome.io.read.nxrsurface(tc.store, s);
                P = double(single(S.Vertices));                              % as the editor receives them, as modes_input does
                tc.verifyEqual(P, py("P", 3), [at ': positions']);
                tc.verifyEqual(S.Faces, py("F", 3) + 1, [at ': faces']);
                tc.verifyEqual(S.Wall, i_read(fullfile(o, 'wall.f64')) ~= 0, [at ': medial wall']);

                [X, G, keep] = rheome.geom.glue(P, S.Faces, S.Wall);
                tc.verifyEqual(G, py("G", 3) + 1, [at ': glued faces']);
                tc.verifyEqual(keep - 1, i_read(fullfile(o, 'keep.f64')), [at ': kept vertices']);
                tc.verifyEqual(X, py("V", 3), [at ': glued positions, bit for bit']);
                tc.verifyEqual(string(rheome.geom.meshkey(G)), string(strtrim(fileread(fullfile(o, 'key.txt')))), [at ': mesh key']);

                H = rheome.operators.hodge_dirac(X, G);
                tc.verifyEqual(H.mesh.F - 1, i_read(fullfile(o, 'iF.f64')), [at ': flips']);
                tc.verifyEqual(H.star1Clamped, 0, [at ': no edge weight floored (the editor floors none)']);

                sol = rheome.eigen.hodge_dirac(H, 200);
                el = i_read(fullfile(o, 'el.f64')); cl = i_read(fullfile(o, 'cl.f64'));
                tc.assertEqual([numel(sol.exactLambda) numel(sol.coexactLambda)], [numel(el) numel(cl)], [at ': the same band']);
                tc.verifyEqual(sol.exactLambda, el, 'RelTol', 1e-9, [at ': exact lambda']);
                tc.verifyEqual(sol.coexactLambda, cl, 'RelTol', 1e-9, [at ': coexact lambda']);
                tc.verifyLessThan(i_modegap(sol.exactModes, py("em", H.mesh.nV)', H.star0), 1e-6, [at ': exact modes']);
                tc.verifyLessThan(i_modegap(sol.coexactModes, py("cm", H.mesh.nF)', H.star2), 1e-6, [at ': coexact modes']);
                clear c
            end
        end

        function chainInputsAreTheEditors(tc)
            % the spectral chain's inputs as the editor reads them: each kernel that carries fused weights, its recording, one 2^15-sample
            % chunk of it; every array against zarr-python's read; and the weights' basis rebuilt here (glue, flips, solve) against the
            % eigenvalues and glued mesh the weights record
            here = fileparts(mfilename('fullpath'));
            for k = 1:numel(tc.subjects)
                s = tc.subjects(k); at = sprintf('subject %d', k);
                [~, out] = system(sprintf('sqlite3 -readonly "%s" "SELECT DISTINCT o.name FROM operator o JOIN subject u ON u.id = o.subject_id JOIN operator w ON json_extract(w.params_json, ''$.kernel'') = o.id WHERE u.name = ''%s'' AND o.kind = ''inverse kernel'' AND w.kind = ''hodge-dirac weights'' ORDER BY o.name"', ...
                    fullfile(tc.store, 'dataset.sqlite'), s));
                kernels = splitlines(strtrim(string(out))); kernels = kernels(kernels ~= "");
                if isempty(kernels), continue, end
                K = rheome.io.read.nxrkernel(tc.store, s, kernels(1));
                tc.assertNotEmpty(K.Recording, [at ': the kernel''s recording']);
                R = rheome.io.read.nxrrecording(tc.store, s, K.Recording);                  % its length first
                R = rheome.io.read.nxrrecording(tc.store, s, K.Recording, Samples=[1 min(2^15, R.nSamples)]);
                tc.verifyTrue(all(K.Channels >= 1 & K.Channels <= numel(R.Channels)) && numel(K.Channels) == size(K.W, 2), [at ': kernel columns are recording channels']);
                tc.verifyGreaterThan(R.Fs, 0);
                tc.verifyEqual(double(R.Data), i_zarr(tc, here, R.Path, size(R.Data), [1 size(R.Data,1) R.Samples]), [at ': a recording chunk']);
                tc.verifyEqual(double(K.W), i_zarr(tc, here, K.Path, size(K.W), []), [at ': the kernel']);
                % THE RESAMPLED COPY (nxr-convert resample: scipy resample_poly of the whole recording): the recording's first samples,
                % 8 channels, resampled here; the copy's ends match from the recording's start, so only this chunk's far end is left out
                r3 = nxrsql300(tc, s, K.Recording);
                if r3 ~= ""
                    R3 = rheome.io.read.nxrrecording(tc.store, s, r3); q = round(R.Fs / R3.Fs);
                    x = rheome.io.read.zarr(R.Path, {[1 8], [1 q * 4096]}); y3 = double(rheome.io.read.zarr(R3.Path, {[1 8], [1 4096]}));
                    y = rheome.filters.resample(double(x), R.Fs, R3.Fs);
                    tc.verifyLessThan(max(abs(y(:, 1:4000) - y3(:, 1:4000)), [], 'all') / max(abs(y3(:, 1:4000)), [], 'all'), 1e-6, [at ': the resampled copy']);
                end
                S = rheome.io.read.nxrsurface(tc.store, s);
                [X, G, keep] = rheome.geom.glue(double(single(S.Vertices)), S.Faces, S.Wall);
                own = K.Weights(strcmp({K.Weights.mesh}, rheome.geom.meshkey(G)));          % the weights fused on this subject's glued mesh
                tc.assertNotEmpty(own, [at ': fused weights on the glued mesh rheome builds']);
                tc.verifyEqual(double(own(1).B), i_zarr(tc, here, own(1).Path, size(own(1).B), []), [at ': the fused weights']);
                Hd = rheome.operators.hodge_dirac(X, G); sol = rheome.eigen.hodge_dirac(Hd, own(1).K + 100);
                mu = sort(sqrt([sol.exactLambda; sol.coexactLambda]));
                tc.assertGreaterThanOrEqual(numel(mu), own(1).K, [at ': the band covers the weights'' modes']);
                tc.verifyEqual(mu(1:own(1).K), own(1).mu, 'RelTol', 1e-9, [at ': the weights'' eigenvalues, rebuilt']);
                % THE FUSION: B rebuilt from the kernel and this basis, against the stored B. A mode's sign is arbitrary (its a and b flip
                % together) and modes whose eigenvalues nearly coincide may mix: the whole matrix, and the typical row, are held
                pkg = rheome.eigen.hodge_package(sol, Hd); Kk = own(1).K;
                pkg.mu = pkg.mu(1:Kk); pkg.kind = pkg.kind(1:Kk); pkg.P = pkg.P(:, 1:Kk);
                B = rheome.inverse.hodge_weights(K.W, pkg, X, G, keep, double(single(S.Vertices)));
                Bs = double(own(1).B); a = 1:Kk; b = Kk + a;
                sg = sign(sum(B(a, :) .* Bs(a, :), 2) + sum(B(b, :) .* Bs(b, :), 2)); B = B .* [sg; sg];
                rel = vecnorm(B - Bs, 2, 2) ./ max(vecnorm(Bs, 2, 2), eps);
                tc.verifyLessThan(norm(B - Bs, 'fro') / norm(Bs, 'fro'), 1e-4, [at ': the fused weights, rebuilt']);
                tc.verifyLessThan(median(rel(b)), 1e-6, [at ': a typical flow row, rebuilt (float32 storage)']);
            end
        end
    end
end

function r = nxrsql300(tc, s, rec)
% the recording resampled from `rec` (its producer names it), '' when none
[~, out] = system(sprintf('sqlite3 -readonly "%s" "SELECT f.name FROM field f JOIN subject u ON u.id = f.subject_id WHERE u.name = ''%s'' AND f.kind = ''recording'' AND json_extract(f.producer_json, ''$.resampled_from.name'') = ''%s'' LIMIT 1"', ...
    fullfile(tc.store, 'dataset.sqlite'), s, rec));
r = string(strtrim(out));
end

function x = i_zarr(tc, here, folder, sz, rng)
% zarr-python's read of the same part
o = [tempname '.f64'];
[st, out] = system(sprintf('uv run --frozen --project "%s" python -I "%s" "%s" "%s" %s 2>&1', tc.nxr, fullfile(here, 'zarr_slice.py'), folder, o, num2str(rng)));
tc.assertEqual(st, 0, out);
fid = fopen(o); x = fread(fid, inf, 'double=>double'); fclose(fid); delete(o);
x = reshape(x, fliplr(sz))';
end

function g = i_modegap(A, B, w)
% the largest mass-norm difference between two sets of unit modes, each sign aligned
s = sign(sum(A .* w .* B, 1)); D = A - B .* s;
g = max(sqrt(sum(D .* w .* D, 1)));
end

function x = i_read(p)
fid = fopen(p); x = fread(fid, inf, 'double=>double'); fclose(fid);
end
