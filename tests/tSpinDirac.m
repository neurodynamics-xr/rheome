classdef tSpinDirac < matlab.unittest.TestCase
% rheome.operators.spin_dirac against the nxr-cortical-flow export (spinDirac.ts) on bunny-closed:
% the three 4F x 4F matrices and H agree entry by entry. Skipped when the reference export is absent.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function matchesCorticalFlowExportOnBunny(tc)
            ws = fullfile(getenv('HOME'), 'workspace');
            d  = fullfile(ws, 'research/results/dirac-eigenmodes-bunny-closed/exports/export-spin');
            fx = fullfile(ws, 'research/code/nxr-cortical-flow/editor/fixtures/bunny-closed');
            tc.assumeTrue(isfile([fx '.faces.u32']) && isfile(fullfile(d, 'bunny.H.f64')));
            V = reshape(i_read([fx '.positions.f32'], 'single'), 3, [])';
            F = reshape(i_read([fx '.faces.u32'], 'uint32'), 3, [])' + 1;
            S = rheome.operators.spin_dirac(V, F);
            for nm = ["extrinsic" "intrinsicHY" "intrinsicYe"]
                p = fullfile(d, "bunny." + nm);
                ip = i_read(p + ".indptr.i32", 'int32'); ix = i_read(p + ".indices.i32", 'int32'); dt = i_read(p + ".data.f64", 'double');
                n = numel(ip) - 1;  R = sparse(repelem((1:n)', diff(ip)), ix + 1, dt, n, n);
                tc.verifyLessThan(norm(S.(nm) - R, 'fro') / norm(R, 'fro'), 1e-12, nm);
                tc.verifyLessThan(norm(S.(nm) - S.(nm)', 'fro'), 1e-12 * norm(R, 'fro'), nm + " symmetric");
            end
            tc.verifyEqual(S.H, i_read(fullfile(d, 'bunny.H.f64'), 'double'), 'AbsTol', 1e-14);
        end

        function vertexToFaceMatchesExport(tc)
            ws = fullfile(getenv('HOME'), 'workspace');
            d  = fullfile(ws, 'research/results/dirac-eigenmodes-bunny-closed/exports/export-spin');
            fx = fullfile(ws, 'research/code/nxr-cortical-flow/editor/fixtures/bunny-closed');
            tc.assumeTrue(isfile([fx '.faces.u32']) && isfile(fullfile(d, 'bunny.Acov.data.f64')));
            V = reshape(i_read([fx '.positions.f32'], 'single'), 3, [])';
            F = reshape(i_read([fx '.faces.u32'], 'uint32'), 3, [])' + 1;
            for nm = ["Acov" "Aplain"]
                p = fullfile(d, "bunny." + nm);
                ip = i_read(p + ".indptr.i32", 'int32'); ix = i_read(p + ".indices.i32", 'int32'); dt = i_read(p + ".data.f64", 'double');
                R = sparse(repelem((1:numel(ip)-1)', diff(ip)), ix + 1, dt, numel(ip)-1, 4*size(V,1));
                A = rheome.operators.vertex_to_face(V, F, nm == "Acov");
                tc.verifyLessThan(norm(A - R, 'fro') / norm(R, 'fro'), 1e-12, nm);
            end
        end

        function extrinsicKillsConstants(tc)
            [V, F] = rheome.geom.icosphere(3);
            S = rheome.operators.spin_dirac(V, F);
            c = repmat([1 0 0 0]', size(F,1), 1);
            tc.verifyLessThan(norm(S.extrinsic*c), 1e-9*norm(S.intrinsicHY*c));
        end
    end
end

function x = i_read(p, t)
fid = fopen(p); x = double(fread(fid, inf, [char(t) '=>double'])); fclose(fid);
end
