% EIGEN  Spectral basis of the surface operators.
%
% Eigenmodes of the Laplace-Beltrami operator, forming the basis every atom is
% projected onto and reconstructed from.
%
% Functions:
%   rheome.eigen.modes       - solve L phi = lambda M phi, M-orthonormal, smallest first (LBO / any pencil)
%   rheome.eigen.dirac_frame - relative-Dirac quaternion eigenbasis
%   rheome.eigen.hodge_dirac - Hodge-Dirac modes split into exact and coexact ladders, cut at one lambda
%                       (nxr-convert solve_hodge; on rheome.operators.hodge_dirac of the glued surface)
%   rheome.eigen.hodge_package - those modes as one vertex potential each, merged by mu (the basis the editor reads)
%   rheome.eigen.lift        - ambient-flat quaternion LIFT of a scalar basis (L (x) I3); used for
%                       Dirac-Connectome = lift of the LB-Connectome basis
%   rheome.eigen.connectome  - low modes of the normalized connectome Laplacian (normalized-adjacency
%                       trick); LB-Connectome modes come from rheome.eigen.modes(A,B,K) directly
%
% Author: Diellor Basha, 2026
