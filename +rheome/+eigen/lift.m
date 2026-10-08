function dbasis = lift(Phis, Lams, Ms)
% EIGEN.LIFT  Lift a scalar eigenbasis to an ambient (flat) quaternion Dirac basis.
%
%   dbasis = rheome.eigen.lift(Phis, Lams, Ms)     % scalar modes + eigenvalues + mass
%   dbasis = rheome.eigen.lift(basis)              % a basis struct with .Phi .Lambda .Mass
%
% Each scalar mode phi_k becomes THREE quaternion modes -- phi_k placed in the imaginary x, y, z
% slots (w = 0), eigenvalue lam_k each -- and the mass is lifted Mq = kron(Ms, I4). This realizes
% L (x) I3 (the ambient, component-wise operator): a flat vector basis in which each cortical
% harmonic points independently along the three ambient axes. Applied to the LB-Connectome scalar
% basis it is the DIRAC-CONNECTOME vector basis (fiber-spread current). Faithful port of
% bst_lift_connectome_dirac.
%
% INPUTS:
%   Phis [n x K] scalar eigenbasis;  Lams [K x 1] eigenvalues;  Ms [n x n] scalar mass
%   (or a single struct with those as .Phi .Lambda .Mass)
%
% OUTPUT (struct dbasis, same layout as rheome.eigen.dirac_frame):
%   .Phi [4n x 3K] quaternion modes ([w;x;y;z] interleaved per vertex)  .Lambda [3K x 1]
%   .Mass kron(Ms,I4)  .nVert n  .nModes 3K
%
% See also: rheome.eigen.modes, rheome.eigen.dirac_frame, rheome.operators.lb_connectome, rheome.filters.tovec
%
% Author: Diellor Basha, 2026 (port of bst_lift_connectome_dirac)

    if nargin == 1 && isstruct(Phis)
        b = Phis;  Lams = b.Lambda;  Ms = b.Mass;  Phis = b.Phi;
    end
    Phis = double(Phis);  [n, K] = size(Phis);
    Phiq = zeros(4*n, 3*K);
    Lamq = zeros(3*K, 1);
    for k = 1:K
        for d = 1:3                                  % d=1,2,3 -> imag x,y,z at rows (v-1)*4 + (d+1)
            c = (k-1)*3 + d;
            Phiq((0:n-1)*4 + (d+1), c) = Phis(:, k);
            Lamq(c) = Lams(k);
        end
    end
    dbasis = struct('Phi', Phiq, 'Lambda', Lamq, 'Mass', kron(Ms, speye(4)), 'nVert', n, 'nModes', 3*K);
end

% Author: Diellor Basha, 2026
