classdef tDiracExportParity < matlab.unittest.TestCase
% rheome.operators.dirac_intrinsic (Crane's intrinsic D) against the nxr-cortical-flow export
% (nxr-compute diracIntrinsicD) on bunny-closed: same sparse 4F x 4V matrix, entry by entry.
% Skipped when the reference export is absent.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function intrinsicDMatchesExport(tc)
            [V, F, d] = i_bunny(tc);
            D = rheome.operators.dirac_intrinsic(V, F);
            R = i_coo(fullfile(d, 'diracIntrinsicD'), size(D));
            tc.verifyLessThan(norm(D - R, 'fro') / norm(R, 'fro'), 1e-12);
        end
        function extrinsicDReport(tc)
            % diracD is the relative D_N with vertex normals; compared through the energy D'*A*D.
            [V, F, d] = i_bunny(tc);
            E = rheome.operators.dirac_extrinsic(V, F, i_angleNormals(V, F));  % cortical-flow's Gauss map
            R = i_coo(fullfile(d, 'diracD'), [4*size(F,1) 4*size(V,1)]);
            [~, MF] = rheome.operators.dirac_intrinsic(V, F);
            ER = R' * MF * R;
            tc.verifyLessThan(norm(E - ER, 'fro') / norm(ER, 'fro'), 1e-6);
        end
    end
end

function [V, F, d] = i_bunny(tc)
ws = fullfile(getenv('HOME'), 'workspace');
d  = fullfile(ws, 'research/results/dirac-eigenmodes-bunny-closed/exports/export');
fx = fullfile(ws, 'research/code/nxr-cortical-flow/editor/fixtures/bunny-closed');
tc.assumeTrue(isfile([fx '.faces.u32']) && isfile(fullfile(d, 'diracD.row.i32')));
V = reshape(i_read([fx '.positions.f32'], 'single'), 3, [])';
F = reshape(i_read([fx '.faces.u32'], 'uint32'), 3, [])' + 1;
end
function R = i_coo(p, sz)
R = sparse(i_read(p + ".row.i32", 'int32') + 1, i_read(p + ".col.i32", 'int32') + 1, i_read(p + ".data.f64", 'double'), sz(1), sz(2));
end
function x = i_read(p, t)
fid = fopen(p); x = double(fread(fid, inf, [char(t) '=>double'])); fclose(fid);
end
function N = i_angleNormals(V, F)
% corner-angle-weighted vertex normals (geometry-central computeVertexNormals), as relativeDirac.ts
N = zeros(size(V));  fn = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);  fn = fn ./ vecnorm(fn, 2, 2);
for k = 1:3
    a = V(F(:,mod(k,3)+1),:) - V(F(:,k),:);  b = V(F(:,mod(k+1,3)+1),:) - V(F(:,k),:);
    ang = atan2(vecnorm(cross(a, b, 2), 2, 2), dot(a, b, 2));
    for c = 1:3, N(:,c) = N(:,c) + accumarray(F(:,k), ang .* fn(:,c), [size(V,1) 1]); end
end
N = N ./ vecnorm(N, 2, 2);
end
