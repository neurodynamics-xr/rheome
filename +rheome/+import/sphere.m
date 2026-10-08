function file = sphere(name, hemi, K, Kc)
% IMPORT.SPHERE  Equip a subject's FreeSurfer registration sphere with every operator + eigenmode.
%
%   file = rheome.import.sphere(name)                 % hemi 'L', K=800 Dirac, Kc=1000 connection
%   file = rheome.import.sphere(name, hemi, K, Kc)
%
% The cortex in +data/<name> is registered (FreeSurfer) to an ico5 sphere of radius 100 mm,
% stored per hemisphere as S.Sphere (Reg.Sphere). This takes that EXACT sphere -- no synthetic
% mesh, so the dimensions are real cortical millimetres and the sphere<->cortex vertex map is
% exact -- and computes on it: the cotan LBO (L,M), the relative-Dirac eigenbasis, the
% Laplace–Beltrami eigenbasis, and the connection-Laplacian operator + eigenmodes. EVERYTHING is
% bundled into ONE struct 'sphere' and cached to +data/<name>/sphere.mat, so a single
% rheome.load.sphere(name) returns the whole sandbox at cortical resolution.
%
% Why the registration sphere and not rheome.geom.icosphere: its LBO spectrum is the clean spherical-
% harmonic sequence (0, sqrt2 x3, sqrt6 x5, ...) to 3 sig figs -- the registration warp is well
% inside FEM tolerance (min triangle angle ~27 deg) -- AND every vertex maps to a real cortical
% location, so a feature found on the sphere maps straight back to the brain.
%
% INPUTS:
%   name  subject cache name (needs rheome.import.surface first, with Reg.Sphere present)
%   hemi  'L' (default) or 'R';  K  Dirac modes (default 800);  Kc  connection modes (default 1000)
%
% OUTPUT struct 'sphere' (saved in sphere.mat):
%   .S       surface struct ON THE SPHERE (.Vertices @ 0.1 m, .Faces, .nV, .VertNormals, .Hemi={})
%   .L .M    LBO stiffness + mass            .basis  Laplace–Beltrami eigenbasis
%   .dirac   relative-Dirac eigenbasis (dbasis)
%   .conn    struct: .C (connection operator A,B,e1,e2,normal,Rt) .Psi .Lambda (eigenmodes)
%   .cortex  [nV x 3] real cortex positions (map results back / visualise on the brain)
%   .globalVertices  index of each sphere vertex into the subject's full cortex
%   .radius .hemi .source .meta   provenance
%
% See also: rheome.import.surface, rheome.operators.connection_laplacian, rheome.eigen.dirac_frame, rheome.load.sphere
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
        error('import:sphere:noSurface', 'No surface cached for ''%s''. Run rheome.import.surface first.', name);
    end
    S = getfield(builtin('load', sf, 'S'), 'S');
    if isempty(S.Sphere)
        error('import:sphere:noReg', 'Surface ''%s'' has no Reg.Sphere (FreeSurfer registration sphere).', name);
    end

    Sh = rheome.utils.hemisphere(S, [], hIdx);              % hemisphere sub-surface: local faces + GlobalVertices
    gv = Sh.GlobalVertices;
    Vs = S.Sphere(gv, :);  Fl = double(Sh.Faces);  nV = size(Vs, 1);
    radius = mean(vecnorm(Vs, 2, 2));
    fprintf('rheome.import.sphere[%s]: %s registration sphere -- %d verts, %d faces, radius %.1f mm\n', ...
        name, hemi, nV, size(Fl,1), 1000*radius);

    Nrm  = Vs ./ vecnorm(Vs, 2, 2);                  % outward sphere normals
    Ssph = struct('Vertices',Vs, 'Faces',Fl, 'nV',nV, 'nF',size(Fl,1), 'VertNormals',Nrm, ...
                  'Hemi',{{}}, 'HemiLabel',{{}}, ...
                  'Comment',sprintf('%s %s reg-sphere (ico5 @ %.0fmm)', name, hemi, 1000*radius), ...
                  'SurfaceFile','', 'Sphere',Vs, 'VertConn',[]);

    [L, M] = rheome.operators.laplace_beltrami(Vs, Fl, 'galerkin');

    fprintf('  Dirac eigenbasis (K=%d) ...\n', K);  t = tic;
    dbasis = rheome.eigen.dirac_frame(Vs, Fl, tau, K, Nrm, {(1:nV)'});
    basis  = rheome.eigen.modes(L, M, K);
    fprintf('    %.1fs\n', toc(t));

    fprintf('  connection Laplacian + eigenmodes (Kc=%d) ...\n', Kc);  t = tic;
    Cc = rheome.operators.connection_laplacian(Vs, Fl);
    [Psi, Lm]   = rheome.eigen.smallest(Cc.A, Cc.B, Kc);   % negative-sigma shift (not naive smallestabs)
    [Lambda, o] = sort(real(diag(Lm)));  Psi = Psi(:, o);
    fprintf('    %.1fs\n', toc(t));

    sphere = struct();
    sphere.S      = Ssph;
    sphere.L      = L;   sphere.M = M;
    sphere.basis  = basis;
    sphere.dirac  = dbasis;
    sphere.conn   = struct('C', Cc, 'Psi', Psi, 'Lambda', Lambda);
    sphere.cortex = Sh.Vertices;                     % real cortex positions (map back)
    sphere.globalVertices = gv;                      % sphere vertex -> subject cortex vertex
    sphere.radius = radius;   sphere.hemi = char(hemi);
    sphere.source = sprintf('%s Reg.Sphere hemi %s', name, hemi);
    sphere.meta   = struct('tau',tau, 'K',K, 'Kc',Kc, 'radius',radius, 'builtFrom',char(name), 'nSub',5);  %#ok<STRNU> saved by name below

    file = fullfile(dsdir, 'sphere.mat');
    builtin('save', file, 'sphere', '-v7.3');
    fprintf('rheome.import.sphere[%s]: bundled operators + eigenmodes -> %s\n', name, file);
end

% Author: Diellor Basha, 2026
