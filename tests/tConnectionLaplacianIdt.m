classdef tConnectionLaplacianIdt < matlab.unittest.TestCase
% rheome.operators.connection_laplacian_idt: on a surface that is already Delaunay it is the same operator
% as rheome.operators.connection_laplacian (spectrum equal, each eigenvalue twice in the real 2n form);
% on a cortex with thousands of flips it stays symmetric and positive semidefinite, and the signposts
% stay inside each vertex's cone.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function sameSpectrumAsInputMeshOperatorOnDelaunaySphere(tc)
            [V, F] = rheome.geom.icosphere(2);
            L = rheome.operators.connection_laplacian_idt(V, F);
            C = rheome.operators.connection_laplacian(V, F);
            tc.verifyLessThan(norm(L - L', 'fro'), 1e-12);
            Cf = full(C.A); ev = sort(real(eig((Cf + Cf')/2)));
            el = sort(eig(full(L)));
            tc.verifyLessThan(max(abs(el(1:2:end) - ev)), 1e-9 * max(1, max(ev)));
            tc.verifyLessThan(max(abs(el(2:2:end) - ev)), 1e-9 * max(1, max(ev)));
        end
        function cortexIsPsdAfterManyFlips(tc)
            d = fullfile(getenv('HOME'), 'workspace/research/results/dirac-eigenmodes-cortex/exports/export-cortex/glued.');
            tc.assumeTrue(isfile([d 'P.f64']));
            V = reshape(i_read([d 'P.f64'], 'double'), 3, [])';
            F = reshape(i_read([d 'F.i32'], 'int32'), 3, [])' + 1;
            [L, M, S] = rheome.operators.connection_laplacian_idt(V, F);
            tc.verifyEqual(S.flips, 10678);
            tc.verifyLessThan(norm(L - L', 'fro'), 1e-9 * norm(L, 'fro'));
            lam = eigs(L, spdiags(M, 0, numel(M), numel(M)), 6, 'smallestabs');
            tc.verifyGreaterThan(min(lam), 0);   % definite: a cortex has no parallel field
            tc.verifyGreaterThanOrEqual(min(S.angle), 0);
            tc.verifyLessThan(max(S.angle ./ S.angleSum(S.F)), 1 + 1e-9);
        end
    end
end
function x = i_read(p, t)
fid = fopen(p); x = double(fread(fid, inf, [char(t) '=>double'])); fclose(fid);
end
