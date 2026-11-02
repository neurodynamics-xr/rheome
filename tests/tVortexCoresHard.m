classdef tVortexCoresHard < matlab.unittest.TestCase
    % rheome.cores against nxr's hard-case goldens (tests/goldens/vortex-hard, synthetic): greedy competition, distance ties,
    % frame gaps, recurrence, reversal and the float32 threshold bin edge. nxr is the reference.
    %
    % Each scenario's curl_heat10 [V x T] goes straight to stats (no modes, no flow), on the first set's mesh (mesh_ref).
    %
    % Author: Diellor Basha, 2026

    properties
        Mesh   % struct: F, V_mm, area_mm2 from the first set
        A
        Dir
        Names
    end

    properties (TestParameter)
        scen = {'competition_single', 'competition_greedy', 'competition_parallel', 'tie', 'gap_zero_frame', ...
                'gap_missing_frames', 'recurrence', 'reversal', 'threshold_edge'}
    end

    methods (TestClassSetup)
        function load(t)
            t.Dir = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'vortex-hard');
            d0 = fullfile(fileparts(mfilename('fullpath')), 'goldens', 'vortex');
            m0 = jsondecode(fileread(fullfile(d0, 'manifest.json')));
            for nm = {'F', 'V_mm', 'area_mm2'}
                t.Mesh.(nm{1}) = tVortexCoresHard.readArray(d0, m0.arrays.(nm{1}));
            end
            t.A = rheome.cores.adjacency(t.Mesh.F + 1, size(t.Mesh.V_mm, 1));
        end
    end

    methods (Static)
        function x = readArray(d, a)
            sh = double(a.shape(:)).';
            f = fopen(fullfile(d, a.file), 'r'); x = fread(f, prod(sh), ['*' strrep(a.dtype, 'float64', 'double')]); fclose(f);
            if numel(sh) == 1, x = x(:); else, x = permute(reshape(x, fliplr(sh)), numel(sh):-1:1); end
        end

        function [G, M, R, steps] = runScenario(t, name)
            d = fullfile(t.Dir, name);
            M = jsondecode(fileread(fullfile(d, 'manifest.json')));
            nm = fieldnames(M.arrays);
            for i = 1:numel(nm), G.(nm{i}) = tVortexCoresHard.readArray(d, M.arrays.(nm{i})); end
            P = jsondecode(fileread(fullfile(t.Dir, 'manifest.json'))).parameters;
            area = t.Mesh.area_mm2;
            if isfield(G, 'area_used_mm2'), area = G.area_used_mm2; end
            [R, steps] = rheome.cores.stats(single(G.curl_heat10), t.A, area, t.Mesh.V_mm, P.f0_hz, ...
                Fs=P.fs_hz, Frames=double(G.frames_used(:)).' + 1);
        end
    end

    methods (Test)
        function allOutputs(t, scen)
            [G, M, R, s] = tVortexCoresHard.runScenario(t, scen);
            p = R.patches;
            % patches: exact; frame is the POSITION in frames_used
            t.verifyEqual(p.frame - 1, double(G.patch_frame), 'patch_frame');
            t.verifyEqual(p.sense, double(G.patch_sense), 'patch_sense');
            t.verifyEqual(p.centre, G.patch_centre_mm, RelTol=1e-6, AbsTol=1e-6);
            t.verifyEqual(p.area, G.patch_area_mm2, RelTol=1e-12);
            % steps: the frame column carries the true frame number
            s(:, 9) = s(:, 9) - 1;
            t.verifyEqual(size(s), size(G.steps), 'steps size');
            t.verifyEqual(s(:, [1 8 9]), G.steps(:, [1 8 9]), 'steps sense/track/frame');
            t.verifyEqual(s(:, 2:7), G.steps(:, 2:7), RelTol=1e-6, AbsTol=1e-6);
            % scalars (null = NaN)
            sc = M.scalars; nm = fieldnames(sc);
            for i = 1:numel(nm)
                if ~isfield(R, nm{i}), continue, end
                v = sc.(nm{i}); if isempty(v), v = NaN; end
                t.verifyEqual(R.(nm{i}), v, 'RelTol', 1e-6, 'AbsTol', 1e-9, nm{i});
            end
            % histograms and occupancy
            t.verifyEqual(double(R.life_hist), double(G.life_hist), 'life_hist');
            t.verifyEqual(double(R.speed_hist), double(G.speed_hist), 'speed_hist');
            t.verifyEqual(double(R.dens_ccw), G.dens_ccw, RelTol=1e-6, AbsTol=1e-7);
            t.verifyEqual(double(R.dens_cw), G.dens_cw, RelTol=1e-6, AbsTol=1e-7);
        end

        function floorThreshold(t, scen)
            % with fewer than 2 % nonzero entries the threshold is the histogram floor lo = 4 max 1e-6 (not threshold_edge)
            if strcmp(scen, 'threshold_edge'), return, end
            [G, ~, R] = tVortexCoresHard.runScenario(t, scen);
            t.verifyEqual(double(R.threshold), 4e-6 * double(max(abs(single(G.curl_heat10(:, double(G.frames_used) + 1))), [], 'all')), RelTol=1e-6);
        end

        function thresholdEdgeFloat32Bin(t)
            % nxr's float32 log-binning: bin 2501, threshold 0.0369453; float64 would give bin 2500, 0.0368208
            [G, M, R] = tVortexCoresHard.runScenario(t, 'threshold_edge');
            sp = M.scenario_parameters;
            C = single(G.curl_heat10);
            [thr, h] = rheome.cores.threshold(C);
            t.verifyEqual(thr, M.scalars.threshold, RelTol=1e-9);
            t.verifyEqual(thr, sp.grid_float32, RelTol=1e-9);
            t.verifyEqual(R.threshold, thr);
            t.verifyGreaterThan(abs(thr - sp.grid_float64) / sp.grid_float64, 3e-3);   % one grid step (0.34 %) above float64
            t.verifyNotEqual(thr, sp.grid_float64);
            t.verifyEqual(find(cumsum(h) / sum(h) >= 0.98, 1) - 1, sp.bin_float32);    % histogram index of the 98th percentile
            % e sits below both thresholds: no core
            t.verifyLessThan(sp.e, sp.grid_float64);
            t.verifyEqual(nnz(R.patches.frame == 12), 0);   % frame 11 (0-based) holds only e
        end
    end
end
