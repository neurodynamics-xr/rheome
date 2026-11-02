classdef tVortexCores < matlab.unittest.TestCase
    % rheome.cores (the threshold-on-curl core detector) against the goldens nxr wrote (tests/goldens/vortex, synthetic).
    %
    % Author: Diellor Basha, 2026

    properties
        G      % struct of golden arrays
        M      % manifest
        R      % rheome.cores.stats result
        steps
        C      % curl movie from the golden curl_heat10 (single)
        A
    end

    methods (TestClassSetup)
        function load(t)
            d = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'vortex');
            t.M = jsondecode(fileread(fullfile(d, 'manifest.json')));
            nm = fieldnames(t.M.arrays);
            for i = 1:numel(nm)
                a = t.M.arrays.(nm{i}); sh = double(a.shape(:)).';
                f = fopen(fullfile(d, a.file), 'r'); x = fread(f, prod(sh), ['*' strrep(a.dtype, 'float64', 'double')]); fclose(f);
                if numel(sh) == 1, x = x(:); else, x = permute(reshape(x, fliplr(sh)), numel(sh):-1:1); end
                t.G.(nm{i}) = x;
            end
            t.C = single(t.G.curl_heat10);
            t.A = rheome.cores.adjacency(t.G.F + 1, size(t.G.V_mm, 1));
            [t.R, t.steps] = rheome.cores.stats(t.C, t.A, t.G.area_mm2, t.G.V_mm, t.M.parameters.f0_hz, Fs=t.M.parameters.fs_hz);
        end
    end

    methods (Test)
        function adjacencyMatches(t)
            A0 = double(t.G.adjacency_indptr); ix = double(t.G.adjacency_indices) + 1;
            B = sparse(repelem((1:numel(A0)-1).', diff(A0)), ix, true, numel(A0)-1, numel(A0)-1);
            t.verifyEqual(t.A, B);
        end

        function modeCoordinates(t)
            [Z, s] = rheome.cores.modecoords(t.G.P, t.G.mass_mm2, t.G.curl_vertex);
            t.verifyEqual(s, t.G.mode_norm_s, RelTol=1e-12);
            t.verifyLessThan(max(abs(Z - t.G.Z_curl), [], 'all') / max(abs(t.G.Z_curl), [], 'all'), 1e-12);
        end

        function heatOperator(t)
            g = rheome.cores.heatgain(t.G.lam_per_mm2, t.M.parameters.core_sigma_mm);
            H = rheome.cores.operator(t.G.P, t.G.mode_norm_s, t.G.Z_curl, g);
            t.verifyClass(H, 'single');
            t.verifyEqual(double(H), t.G.curl_heat10, 'RelTol', 1e-6, 'AbsTol', 1e-6 * max(abs(t.G.curl_heat10), [], 'all'));
        end

        function hatOperatorLevels(t)
            for l = t.M.parameters.levels(:).'
                g = rheome.cores.hatgain(t.G.lam_per_mm2, sum(t.G.area_mm2), l);
                H = rheome.cores.operator(t.G.P, t.G.mode_norm_s, t.G.Z_curl, g);
                ref = t.G.(sprintf('curl_hat_level%d', l));
                t.verifyEqual(double(H), ref, 'RelTol', 1e-6, 'AbsTol', 1e-6 * max(abs(ref), [], 'all'));
            end
        end

        function thresholdMatches(t)
            [thr, h] = rheome.cores.threshold(t.C);
            t.verifyEqual(thr, t.M.scalars.threshold, RelTol=1e-9);
            t.verifyEqual(thr, t.M.scalars.threshold_loop_reference, RelTol=1e-9);
            t.verifyEqual(sum(h), numel(t.C));
        end

        function patchesExact(t)
            p = t.R.patches;
            t.verifyEqual(p.frame - 1, double(t.G.patch_frame));
            t.verifyEqual(p.sense, double(t.G.patch_sense));
            t.verifyEqual(p.centre, t.G.patch_centre_mm, RelTol=1e-6, AbsTol=1e-6);
            t.verifyEqual(p.area, t.G.patch_area_mm2, RelTol=1e-12);
        end

        function occupancy(t)
            t.verifyEqual(double(t.R.dens_ccw), t.G.dens_ccw, RelTol=1e-6, AbsTol=1e-7);
            t.verifyEqual(double(t.R.dens_cw), t.G.dens_cw, RelTol=1e-6, AbsTol=1e-7);
        end

        function stepsMatch(t)
            ref = t.G.steps; s = t.steps; s(:, 9) = s(:, 9) - 1;
            t.verifyEqual(size(s), size(ref));
            t.verifyEqual(s(:, [1 8 9]), ref(:, [1 8 9]));
            t.verifyEqual(s(:, 2:7), ref(:, 2:7), RelTol=1e-6, AbsTol=1e-6);
        end

        function trackMeasures(t)
            sc = t.M.scalars; nm = fieldnames(sc);
            for i = 1:numel(nm)
                if ~isfield(t.R, nm{i}), continue, end
                t.verifyEqual(t.R.(nm{i}), sc.(nm{i}), 'RelTol', 1e-6, 'AbsTol', 1e-9, nm{i});
            end
        end

        function histograms(t)
            t.verifyEqual(double(t.R.life_hist), double(t.G.life_hist));
            t.verifyEqual(double(t.R.speed_hist), double(t.G.speed_hist));
        end

        function plantedVorticesFound(t)
            % the detector reads back the planted answer: two cores a frame, one each sense
            t.verifyEqual(t.R.cores_per_frame, 2);
            t.verifyEqual(sort(unique(t.R.patches.sense)).', [-1 1]);
        end

        function faceCurlVsNxr(t)
            % rheome.differential.curl on the golden flow against nxr's curl_vertex. nxr's curl is per metre (the golden
            % mesh is stored in mm), so the mesh goes in as metres; nxr first projects the flow onto the tangent plane
            % (vertex normals) and casts positions to float32. Rheome's curl has no projection: reported both ways.
            S = struct('Vertices', t.G.V_mm / 1000, 'Faces', t.G.F + 1);
            nV = size(t.G.V_mm, 1); nT = size(t.G.flow, 3); nrm = t.G.normals;
            Jt = permute(t.G.flow, [1 2 3]); Jp = Jt - sum(Jt .* nrm, 2) .* nrm;
            sc = max(abs(t.G.curl_vertex), [], 'all');
            c0 = rheome.differential.curl(reshape(permute(Jt, [2 1 3]), 3*nV, nT), S);
            c1 = rheome.differential.curl(reshape(permute(Jp, [2 1 3]), 3*nV, nT), S);
            r0 = max(abs(c0 - t.G.curl_vertex), [], 'all') / sc; r1 = max(abs(c1 - t.G.curl_vertex), [], 'all') / sc;
            fprintf('face curl vs nxr curl_vertex (max diff / max |curl|): raw flow %.3g, tangent-projected flow %.3g\n', r0, r1);
            t.verifyLessThan(r1, 1e-6);
            t.verifyLessThan(r0, 1e-3);
        end
    end
end
