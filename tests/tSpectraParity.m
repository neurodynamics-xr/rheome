classdef tSpectraParity < matlab.unittest.TestCase
% Eigenvalue parity with the nxr-cortical-flow Python reference (solve-diracs.py) on bunny-closed:
% the 128 Hodge-Dirac eigenvalues nearest sigma = 1 (pencil K x = mu M x), and the 128 smallest mu of
% Crane's intrinsic D and the relative D_N (squares D' W_F D u = mu^2 M_V u, shift-invert at -1).
% Skipped when the reference results are absent. Small mesh only: a few eigensolves, no large runs.
%
% Author: Diellor Basha, 2026

    properties (Constant)
        Res = fullfile(getenv('HOME'), 'workspace/research/results/dirac-eigenmodes-bunny-closed')
        Fx  = fullfile(getenv('HOME'), 'workspace/research/code/nxr-cortical-flow/editor/fixtures/bunny-closed')
    end
    methods (Test)
        function hodgeSpectrum(tc)
            [V, F] = i_mesh(tc);
            H = rheome.operators.hodge_dirac(V, F);
            mu = eigs(H.K, spdiags(H.M, 0, H.N, H.N), 128, 1.0);
            ref = i_npy(fullfile(tc.Res, 'hodge_mu.npy'));
            tc.verifyLessThan(max(abs(sort(mu) - ref)) / max(abs(ref)), 1e-8);  % absolute: kernel eigenvalues are ~0
        end
        function intrinsicDiracSpectrum(tc)
            [V, F] = i_mesh(tc);
            tc.verifyLessThan(i_sq(tc, rheome.operators.dirac_intrinsic(V, F), V, F, 'diracIntrinsicD_mu.npy'), 1e-6);
        end
    end
end

function [V, F] = i_mesh(tc)
fx = tc.Fx; tc.assumeTrue(isfile([fx '.faces.u32']) && isfile(fullfile(tc.Res, 'hodge_mu.npy')));
V = reshape(i_read([fx '.positions.f32'], 'single'), 3, [])';
F = reshape(i_read([fx '.faces.u32'], 'uint32'), 3, [])' + 1;
end
function err = i_sq(tc, D, V, F, refname)
nF = size(F,1); nV = size(V,1);
e1 = V(F(:,2),:)-V(F(:,1),:); e2 = V(F(:,3),:)-V(F(:,1),:); A = vecnorm(cross(e1,e2,2),2,2)/2;
W = spdiags(repelem(A,4), 0, 4*nF, 4*nF);
MV = spdiags(repelem(accumarray(F(:), repmat(A/3,3,1), [nV 1]), 4), 0, 4*nV, 4*nV);
lam = sort(eigs(D'*W*D, MV, 128, -1.0));
ref = i_npy(fullfile(tc.Res, refname));
err = max(abs(sqrt(max(lam,0)) - ref)) / max(ref);
end
function x = i_npy(p)
fid = fopen(p); fread(fid, 8); hl = fread(fid, 1, 'uint16'); fread(fid, hl); x = fread(fid, inf, 'double'); fclose(fid);
end
function x = i_read(p, t)
fid = fopen(p); x = double(fread(fid, inf, [char(t) '=>double'])); fclose(fid);
end
