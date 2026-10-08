function dbasis = dirac_frame(V, F, tau, nModes, N, hemis, normalize)
% EIGEN.DIRAC_FRAME  Per-hemisphere eigenbasis of the relative Dirac.
%
%   dbasis = rheome.eigen.dirac_frame(V, F, tau, nModes)
%   dbasis = rheome.eigen.dirac_frame(V, F, tau, nModes, N)          % supply the Gauss map
%   dbasis = rheome.eigen.dirac_frame(V, F, tau, nModes, N, hemis)   % supply the hemisphere split
%   dbasis = rheome.eigen.dirac_frame(..., hemis, normalize)         % normalize=false -> physical eigenvalues
%
% Reproduces Brainstorm bst_dirac's basis: split the cortex into hemispheres, build the
% relative Dirac (default: the (1-tau)*intrinsic + tau*extrinsic blend)
% A_h = (1-tau)(D_int^2/sL)+tau(E/sE) per hemisphere
% (rheome.operators.dirac_frame), solve  A_h phi = lambda B_h phi  for nModes each (rheome.eigen.modes;
% M-orthonormalization is the Rayleigh-Ritz step for the quaternion degeneracy), and stack
% the per-hemisphere modes -- matching bst_dirac's [L~_L, L~_R]. Returns the module's Dirac
% dbasis shape, so rheome.filters.impulse / rheome.forward.reconstruct consume it unchanged. This is THE
% Dirac basis of the standalone module -- the proper curvature-aware operator, matching
% Brainstorm; there is no flat-Dirac alternative.
%
% HEMISPHERE SPLIT: pass 'hemis' = {leftVerts, rightVerts} (rheome.io.read.surface .Hemi, from the
% Structures atlas Cortex L/R scouts) -- the AUTHORITATIVE split. Only when it is not given
% does this fall back to connected components (fine for synthetic meshes, fragile on real
% cortex with non-manifold edges).
%
% INPUTS:  V,F surface;  tau in [0,1] (default caller passes 0.5, as Brainstorm);
%          nModes modes PER hemisphere;  N (optional) unit vertex normals;
%          hemis (optional) cell of per-hemisphere vertex-index lists.
% OUTPUT (struct dbasis):
%   .Phi [4V x M] B-orthonormal (M = nModes summed over hemispheres), .Lambda [M x 1],
%   .Mass = kron(Mass_full, I4), .nVert, .nModes = M, .Tau, .Hemisphere [M x 1],
%   .Normalize (the scale convention), .Scales [nComp x 2] = [sL sE] per hemisphere. To read
%   the spectrum in PHYSICAL units without rebuilding:  Lambda_phys = Lambda .* Scales(Hemisphere,1).
%
% See also: rheome.io.read.surface, rheome.operators.dirac_frame, rheome.operators.dirac_intrinsic_sq, bst_dirac
%
% Author: Diellor Basha, 2026

    if nargin < 5, N = []; end
    if nargin < 6, hemis = []; end
    if nargin < 7 || isempty(normalize), normalize = true; end
    nV = size(V,1);
    if ~isempty(hemis)
        comps = cellfun(@(v) v(:), hemis(:)', 'uni', 0);  % authoritative partition (Structures atlas)
    else
        comps = i_components(V, F);                        % fallback: connected components
    end

    Lam = []; Hemi = []; Ii = []; Jj = []; Vv = [];  col0 = 0;
    Sc = zeros(numel(comps), 2);                          % [sL sE] per component (the co-normalization scales)
    for c = 1:numel(comps)
        vidx = comps{c}(:);                               % global vertex indices of this component
        [Vc, Fc, Nc] = i_submesh(V, F, N, vidx);
        [Ac, Bc, sc] = rheome.operators.dirac_frame(Vc, Fc, tau, Nc, normalize);
        Sc(c,:) = sc;                                     % physical lambda = lambda_norm * sL(hemisphere)
        bc = rheome.eigen.modes(Ac, Bc, nModes);                 % [4nVc x mc] B-orthonormal, ascending
        mc = size(bc.Phi, 2);
        % scatter Phi_c (rows over the component's 4-vertex blocks) into the full [4V x .] Phi
        rowGlob = reshape((vidx' - 1)*4 + (1:4)', [], 1); % local row r -> global row
        [ri, ci, vi] = find(bc.Phi);
        Ii = [Ii; rowGlob(ri)];  Jj = [Jj; col0 + ci];  Vv = [Vv; vi];
        Lam  = [Lam; bc.Lambda(:)];
        Hemi = [Hemi; c*ones(mc,1)];
        col0 = col0 + mc;
    end
    Phi = sparse(Ii, Jj, Vv, 4*nV, col0);

    M = rheome.operators.mass(V, F, 'galerkin');
    dbasis = struct('Phi', Phi, 'Lambda', Lam, 'Mass', kron(M, speye(4)), ...
                    'nVert', nV, 'nModes', col0, 'Tau', tau, 'Hemisphere', Hemi, ...
                    'Normalize', normalize, 'Scales', Sc);   % Scales [nComp x 2] = [sL sE]; physical lambda = Lambda.*Sc(Hemisphere,1)
end

% ---- connected components of the mesh graph (hemispheres are the components) ----
% Vertex grouping only (not a geometry-modifying surface split); a clean cortex has
% two disconnected hemispheres, so this returns exactly the lh/rh partition.
function comps = i_components(V, F)
    nV = size(V,1);
    E  = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    G  = graph(E(:,1), E(:,2), [], nV);
    ci = conncomp(G)';
    comps = arrayfun(@(k) find(ci == k), 1:max(ci), 'uni', 0);
end

% ---- extract a component as its own mesh (re-indexed faces) ----
function [Vc, Fc, Nc] = i_submesh(V, F, N, vidx)
    remap = zeros(size(V,1),1);  remap(vidx) = 1:numel(vidx);
    keep = all(ismember(F, vidx), 2);
    Fc = remap(F(keep,:));
    Vc = V(vidx,:);
    if isempty(N), Nc = []; else, Nc = N(vidx,:); end
end

% Author: Diellor Basha, 2026
