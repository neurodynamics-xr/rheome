function file = cortex(name, hemi, K, Kc)
% IMPORT.CORTEX  Equip a subject's FOLDED cortical hemisphere with every operator + eigenmode.
%
%   file = rheome.import.cortex(name)                 % hemi 'L', K=800 Dirac, Kc=1000 connection
%   file = rheome.import.cortex(name, hemi, K, Kc)
%
% The folded-metric TWIN of rheome.import.sphere. Same mesh, real cortical embedding: it takes one
% hemisphere of the subject's cortex (the actual folded 3-D vertex positions) and computes on
% it the cotan LBO (L,M), the relative-Dirac eigenbasis, the Laplace–Beltrami eigenbasis, and
% the connection-Laplacian operator + eigenmodes -- exactly the operators sphere.mat holds, but
% on the FOLDED metric instead of the 100 mm sphere. Bundled into one struct 'cortex' and cached
% to +data/<name>/cortex.mat, so rheome.load.cortex(name) returns the whole thing instantly.
%
% Because the sphere and cortex share connectivity (identical faces, identity vertex map via
% rheome.load.sphere.globalVertices), the two bundles are the SAME mesh in two embeddings -- the basis
% for matched-vs-mismatched filtering (the sphere-vs-cortex comparison, since removed). Named 'cortex' (not 'eigen') to
% avoid shadowing the rheome.eigen.* package.
%
% INPUTS:
%   name  subject cache name (needs rheome.import.surface first);  hemi 'L' (default) or 'R'
%   K     Dirac modes (default 800);  Kc  connection modes (default 1000)
%
% OUTPUT struct 'cortex' (saved in cortex.mat):
%   .S       surface struct on the FOLDED cortex (.Vertices in m, .Faces, .nV, .VertNormals, .Hemi={})
%   .L .M    LBO stiffness + mass            .basis  Laplace–Beltrami eigenbasis
%   .dirac   relative-Dirac eigenbasis (dbasis)
%   .conn    struct: .C (connection operator A,B,e1,e2,normal,Rt) .Psi .Lambda (eigenmodes)
%   .globalVertices  index of each vertex into the subject's full cortex (== sphere's, same order)
%   .hemi .source .meta   provenance
%
% See also: rheome.import.sphere, rheome.operators.connection_laplacian, rheome.eigen.dirac_frame, rheome.load.cortex
%
% Author: Diellor Basha, 2026

    if nargin<2 || isempty(hemi), hemi = 'L';  end
    if nargin<3 || isempty(K),    K    = 800;  end
    if nargin<4 || isempty(Kc),   Kc   = 1000; end
    tau  = 0.5;
    hIdx = 1 + strcmpi(hemi, 'R');

    dsdir = fullfile(rheome.load.root(), char(name));
    sf = fullfile(dsdir, 'surface.mat');
    if ~exist(sf, 'file')
        error('import:cortex:noSurface', 'No surface cached for ''%s''. Run rheome.import.surface first.', name);
    end
    S = getfield(builtin('load', sf, 'S'), 'S');

    Sh = rheome.utils.hemisphere(S, [], hIdx);              % folded hemisphere: local faces + GlobalVertices
    gv = Sh.GlobalVertices;
    Vc = Sh.Vertices;  Fc = double(Sh.Faces);  nV = size(Vc, 1);
    Nrm = Sh.VertNormals;                            % real cortical normals (Dirac extrinsic term)
    if isempty(Nrm) || size(Nrm,1) ~= nV, Nrm = i_vnormals(Vc, Fc); end
    fprintf('rheome.import.cortex[%s]: %s folded cortex -- %d verts, %d faces\n', name, hemi, nV, size(Fc,1));

    Scx = struct('Vertices',Vc, 'Faces',Fc, 'nV',nV, 'nF',size(Fc,1), 'VertNormals',Nrm, ...
                 'Hemi',{{}}, 'HemiLabel',{{}}, 'Comment',sprintf('%s %s cortex (folded)', name, hemi), ...
                 'SurfaceFile','', 'Sphere',[], 'VertConn',[]);

    [L, M] = rheome.operators.laplace_beltrami(Vc, Fc, 'galerkin');

    fprintf('  Dirac eigenbasis (K=%d) ...\n', K);  t = tic;
    dbasis = rheome.eigen.dirac_frame(Vc, Fc, tau, K, Nrm, {(1:nV)'});
    basis  = rheome.eigen.modes(L, M, K);
    fprintf('    %.1fs\n', toc(t));

    fprintf('  connection Laplacian + eigenmodes (Kc=%d) ...\n', Kc);  t = tic;
    Cc = rheome.operators.connection_laplacian(Vc, Fc);
    [Psi, Lm]   = rheome.eigen.smallest(Cc.A, Cc.B, Kc);   % negative-sigma shift (not naive smallestabs)
    [Lambda, o] = sort(real(diag(Lm)));  Psi = Psi(:, o);
    fprintf('    %.1fs\n', toc(t));

    cortex = struct();
    cortex.S      = Scx;
    cortex.L      = L;   cortex.M = M;
    cortex.basis  = basis;
    cortex.dirac  = dbasis;
    cortex.conn   = struct('C', Cc, 'Psi', Psi, 'Lambda', Lambda);
    cortex.globalVertices = gv;                      % vertex -> subject cortex vertex (same order as sphere)
    cortex.hemi   = char(hemi);
    cortex.source = sprintf('%s cortex hemi %s (folded)', name, hemi);
    cortex.meta   = struct('tau',tau, 'K',K, 'Kc',Kc, 'builtFrom',char(name));  %#ok<STRNU> saved by name below

    file = fullfile(dsdir, 'cortex.mat');
    builtin('save', file, 'cortex', '-v7.3');
    fprintf('rheome.import.cortex[%s]: bundled operators + eigenmodes -> %s\n', name, file);
end

% ----- fallback vertex normals (area-weighted), if the surface carries none -----
function N = i_vnormals(V, F)
    fn = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:));   % area-weighted face normals
    N  = [accumarray(F(:), repmat(fn(:,1),3,1), [size(V,1) 1]), ...
          accumarray(F(:), repmat(fn(:,2),3,1), [size(V,1) 1]), ...
          accumarray(F(:), repmat(fn(:,3),3,1), [size(V,1) 1])];
    N  = N ./ max(vecnorm(N,2,2), eps);
end

% Author: Diellor Basha, 2026
