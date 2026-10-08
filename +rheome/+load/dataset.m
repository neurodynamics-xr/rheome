function D = dataset(name, tau, K)
% LOAD.DATASET  Load everything cached for a dataset in one struct.
%
%   D = rheome.load.dataset(name)              % surface (+operators), study & dirac if present (tau=0.5,K=400)
%   D = rheome.load.dataset(name, tau, K)
%
% Returns a struct with whatever is cached:
%   .S .L .M .cortexFile   (always -- surface + operators)
%   .dbasis .basis         (if the (tau,K) Dirac eigenbasis is cached; basis = Laplace–Beltrami)
%   .chan .hm .rec .ncov .studyDir .dataName   (if a study is cached)
%   .inverse               (if a Brainstorm inverse kernel is cached; see rheome.import.inverse)
%   .sphere                (the registration-sphere sandbox: operators + eigenmodes; see rheome.import.sphere)
%   .cortex                (the folded-cortex sandbox: operators + eigenmodes; see rheome.import.cortex)
%   .bases                 (per-hemisphere operators + eigenbases, all families; see rheome.import.bases)
%   .connectome            (whole-brain connectome operators + eigenbases; see rheome.import.connectome)
%   .atlas                 (the surface's parcellations -- the ROI axis; see rheome.import.atlas)
%
% Geometry-only datasets simply have no study fields. This is the convenience entry
% point for demos: import once, then  D = rheome.load.dataset(name)  and go.
%
% See also: rheome.import.dataset, rheome.load.surface, rheome.load.dirac, rheome.load.study
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(tau), tau = 0.5; end
    if nargin < 3 || isempty(K),   K   = 400; end

    sf = fullfile(rheome.load.root(), char(name), 'surface.mat');
    if ~exist(sf, 'file')
        error('load:dataset:missing', ...
            'No dataset ''%s'' cached. Run  rheome.import.dataset(''%s'', studyDir, cortexFile, dataName)  first.', name, name);
    end
    Csurf = builtin('load', sf);
    D.S = Csurf.S;  D.L = Csurf.L;  D.M = Csurf.M;  D.cortexFile = Csurf.cortexFile;

    dfile = fullfile(rheome.load.root(), char(name), sprintf('dirac__tau%.2f__K%d.mat', tau, K));
    if ~exist(dfile, 'file')                        % fall back to any cached K (e.g. the ico sphere's K=800)
        g = dir(fullfile(rheome.load.root(), char(name), sprintf('dirac__tau%.2f__K*.mat', tau)));
        if ~isempty(g), dfile = fullfile(g(1).folder, g(1).name); end
    end
    if exist(dfile, 'file')
        Cd = builtin('load', dfile, 'dbasis', 'basis');  D.dbasis = Cd.dbasis;  D.basis = Cd.basis;
    end

    sfile = fullfile(rheome.load.root(), char(name), 'study.mat');
    if exist(sfile, 'file')
        st = rheome.load.study(name);
        D.chan = st.chan;  D.hm = st.hm;  D.rec = st.rec;  D.ncov = st.ncov;
        D.studyDir = st.studyDir;  D.dataName = st.dataName;
    end

    ifile = fullfile(rheome.load.root(), char(name), 'inverse_bst.mat');
    if exist(ifile, 'file'), D.inverse = rheome.load.inverse(name); end

    if exist(fullfile(rheome.load.root(), char(name), 'sphere.mat'), 'file')
        D.sphere = rheome.load.sphere(name);
    end
    if exist(fullfile(rheome.load.root(), char(name), 'cortex.mat'), 'file')
        D.cortex = rheome.load.cortex(name);
    end
    if exist(fullfile(rheome.load.root(), char(name), 'bases.mat'), 'file')
        D.bases = rheome.load.bases(name, tau, K);        % per-hemisphere operators + eigenbases (all families)
    end
    if exist(fullfile(rheome.load.root(), char(name), 'connectome.mat'), 'file')
        D.connectome = rheome.load.connectome(name);      % whole-brain connectome operators + eigenbases
    end
    if exist(fullfile(rheome.load.root(), char(name), 'atlas.mat'), 'file')
        D.atlas = rheome.load.atlas(name);                % parcellations -- the ROI axis
    end
end

% Author: Diellor Basha, 2026
