classdef tHodgeDirac < matlab.unittest.TestCase
% rheome.operators.intrinsic_delaunay / hodge_dirac against the nxr-cortical-flow export on bunny-closed
% (hodge.{indptr,indices,data} = K, hodge.mass = M), plus d1*d0 = 0 and the Delaunay property.
% Skipped when the reference export is absent.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function matchesExport(tc)
            ws = fullfile(getenv('HOME'), 'workspace');
            d  = fullfile(ws, 'research/results/dirac-eigenmodes-bunny-closed/exports/export');
            fx = fullfile(ws, 'research/code/nxr-cortical-flow/editor/fixtures/bunny-closed');
            tc.assumeTrue(isfile([fx '.faces.u32']) && isfile(fullfile(d, 'hodge.mass.f64')));
            V = reshape(i_read([fx '.positions.f32'], 'single'), 3, [])';
            F = reshape(i_read([fx '.faces.u32'], 'uint32'), 3, [])' + 1;
            H = rheome.operators.hodge_dirac(V, F);
            ip = i_read(fullfile(d,'hodge.indptr.i32'),'int32'); ix = i_read(fullfile(d,'hodge.indices.i32'),'int32');
            dt = i_read(fullfile(d,'hodge.data.f64'),'double'); n = numel(ip)-1;
            R = sparse(repelem((1:n)', diff(ip)), ix+1, dt, n, n);
            tc.verifyEqual(H.N, n);
            tc.verifyLessThan(norm(H.K - R,'fro')/norm(R,'fro'), 1e-9, 'K');
            Mr = i_read(fullfile(d,'hodge.mass.f64'),'double');
            tc.verifyLessThan(norm(H.M - Mr, inf)/norm(Mr, inf), 1e-9, 'M');
        end
        function matchesExportOnGluedCortex(tc)
            % 10 678 flips on a 39 380-face cortex (cortex_pial_low, hemispheres glued): the flip logic proper
            d = fullfile(getenv('HOME'), 'workspace/research/results/dirac-eigenmodes-cortex/exports/export-cortex/glued.');
            tc.assumeTrue(isfile([d 'hodge.mass.f64']));
            V = reshape(i_read([d 'P.f64'], 'double'), 3, [])';
            F = reshape(i_read([d 'F.i32'], 'int32'), 3, [])' + 1;
            H = rheome.operators.hodge_dirac(V, F);
            ip = i_read([d 'hodge.indptr.i32'],'int32'); ix = i_read([d 'hodge.indices.i32'],'int32');
            n = numel(ip)-1;
            R = sparse(repelem((1:n)', diff(ip)), ix+1, i_read([d 'hodge.data.f64'],'double'), n, n);
            tc.verifyEqual(H.mesh.flips, 10678);
            tc.verifyLessThan(norm(H.K - R,'fro')/norm(R,'fro'), 1e-9, 'K');
            Mr = i_read([d 'hodge.mass.f64'],'double');
            tc.verifyLessThan(norm(H.M - Mr, inf)/norm(Mr, inf), 1e-9, 'M');
        end
        function complexAndDelaunay(tc)
            [V, F] = rheome.geom.icosphere(3);
            V = V .* [1 1.3 0.6];                          % an ellipsoid: some edges need flipping
            H = rheome.operators.hodge_dirac(V, F);
            tc.verifyEqual(full(max(abs(H.d1*H.d0), [], 'all')), 0);
            tc.verifyEqual(H.mesh.nonDelaunay, 0);
            tc.verifyGreaterThanOrEqual(min(H.star1), -1e-6);
            tc.verifyLessThan(norm(H.mesh.L*ones(H.mesh.nV,1)), 1e-9);
            tc.verifyEqual(sum(H.mesh.mass), sum(H.star0));
        end
    end
end
function x = i_read(p, t)
fid = fopen(p); x = double(fread(fid, inf, [char(t) '=>double'])); fclose(fid);
end
