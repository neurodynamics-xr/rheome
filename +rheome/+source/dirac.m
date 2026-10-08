function out = dirac(studyDir, cortexFile, dataName, K, measure)
% SOURCE.DIRAC  End-to-end Dirac-eigenmode MEG source mapping (no plotting).
%
%   out = rheome.source.dirac(studyDir, cortexFile, dataName)      % from Brainstorm files
%   out = rheome.source.dirac(studyDir, cortexFile, dataName, K)
%   out = rheome.source.dirac(datasetName)                         % from a cached dataset (+data)
%   out = rheome.source.dirac(datasetName, [], [], K)
%   out = rheome.source.dirac(..., K, measure)                     % InverseMeasure (see below)
%
% Builds the relative-Dirac eigenbasis (rheome.eigen.dirac_frame, tau=0.5, as
% Brainstorm's bst_dirac), projects the leadfield, solves the whitened minimum-norm
% inverse, and returns the result -- including the mode-coefficient time series that
% later dynamics analysis consumes. Faithful to bst_dirac + bst_inverse_dirac. This is
% the compute core of the Dirac route (rheome.flow.context, Method=dirac).
%
% INPUTS (two forms):
%   FILE form:    studyDir (channel_*/noisecov_full/headmodel_surf_os_meg*), cortexFile,
%                 dataName -- read + computed fresh.
%   DATASET form: the first argument is a cached dataset NAME (see rheome.import.dataset); the
%                 surface, operators, Dirac + scalar eigenbases and study are LOADED from
%                 +data (no recompute). Detected via rheome.load.has(studyDir).
%   K : Laplace–Beltrami modes / Dirac modes-per-hemisphere (default 400).
%   measure : InverseMeasure for rheome.inverse.dirac -- 'dspm2018' (default, unitless noise-
%             normalized statistic) | 'amplitude' (whitened min-norm CURRENT, physical
%             units A.m -- use this to read/visualize true dipole magnitudes) | 'sloreta'.
%             Note the mode-coefficient series out.c is measure-INDEPENDENT (STAGE 3);
%             only out.Results.ImagingKernel / out.Jframe depend on it.
%
% OUTPUT (struct out):
%   .S surface   .basis Laplace–Beltrami basis   .dbasis relative-Dirac eigenbasis
%   .Results rheome.inverse.dirac output   .c [nModes x nT] mode-coeff series   .F selected data
%   .Time .sfreq .iSel .SurfaceFile   .frame peak-GFP sample   .Jframe [3nV x 1]
%   .studyDir .dataName   (where the Brainstorm study lives -- for rheome.import.inverse / comparison)
%
% See also: rheome.forward.dirac, rheome.inverse.dirac, rheome.import.dataset
%
% Author: Diellor Basha, 2026

    if nargin < 4 || isempty(K), K = 400; end
    if nargin < 5 || isempty(measure), measure = 'dspm2018'; end

    if nargin >= 1 && rheome.load.has(studyDir)
        % ---- DATASET form: load everything from the +data cache ----
        name = studyDir;
        D = rheome.load.dataset(name, 0.5, K);
        if ~isfield(D, 'dbasis')
            error('source:dirac:noDirac', 'Dataset ''%s'' has no Dirac (tau=0.5, K=%d). Run rheome.import.dirac(''%s'',0.5,%d).', name, K, name, K);
        end
        if ~isfield(D, 'chan')
            error('source:dirac:noStudy', 'Dataset ''%s'' has no study cached. Run rheome.import.study(''%s'', studyDir, dataName).', name, name);
        end
        S = D.S;  basis = D.basis;  dbasis = D.dbasis;
        chan = D.chan;  rec = D.rec;  ncov = D.ncov;  hm = D.hm;
        cortexFile = D.cortexFile;
        studyDirOut = D.studyDir;  dataNameOut = D.dataName;   % where the Brainstorm study lives
    else
        % ---- FILE form: read Brainstorm files + compute (glob the leadfield, e.g. _02) ----
        chanF = i_find(studyDir, 'channel_*.mat');
        ncovF = fullfile(studyDir, 'noisecov_full.mat');
        hmF   = i_find(studyDir, 'headmodel_surf_os_meg*.mat');
        dataF = fullfile(studyDir, dataName);
        for f = {chanF, ncovF, hmF, dataF}
            if isempty(f{1}) || ~exist(f{1},'file'), error('source:dirac:missing','Missing file: %s', f{1}); end
        end
        chan = rheome.io.read.channel(chanF);
        rec  = rheome.io.read.recording(dataF);
        ncov = rheome.io.read.noisecov(ncovF);
        hm   = rheome.io.read.headmodel(hmF);
        S    = rheome.io.read.surface(cortexFile);
        [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces, 'galerkin');
        basis  = rheome.eigen.modes(L, M, K);                    % Laplace–Beltrami (returned for reference)
        dbasis = rheome.eigen.dirac_frame(S.Vertices, S.Faces, 0.5, K, S.VertNormals, S.Hemi);   % atlas L/R split
        studyDirOut = studyDir;  dataNameOut = dataName;
    end

    if hm.nV ~= S.nV
        error('source:dirac:grid', 'Leadfield grid has %d sources but the cortex has %d vertices.', hm.nV, S.nV);
    end

    % --- select good MEG channels (consistent across leadfield / data / noisecov) ---
    isMEG = strcmpi(chan.Type, 'MEG');
    iSel  = find(isMEG(:) & (rec.ChannelFlag(:) == 1));
    Gain  = hm.Gain(iSel, :);
    F     = []; if ~isempty(rec.F), F = rec.F(iSel, :); end   % empty for raw/kernel-only studies
    Cnoise = ncov.NoiseCov(iSel, iSel);
    chTypes = chan.Type(iSel);

    % --- forward + inverse (data-independent: the kernel does NOT need the recording) ---
    Gm  = rheome.forward.dirac(Gain, dbasis);
    Res = rheome.inverse.dirac(Gm, dbasis, struct('NoiseCov', Cnoise, ...
              'FourthMoment', i_sub(ncov.FourthMoment, iSel), 'nSamples', i_sub(ncov.nSamples, iSel)), ...
              struct('ChannelTypes', {chTypes}, 'InverseMeasure', measure));

    % --- apply: mode coefficients + a peak-GFP display frame (only when data is present) ---
    if isempty(F)
        c = [];  frame = [];  Jframe = [];              % kernel-only: streaming happens downstream
    else
        c   = Res.ImagingKernelMode * F;                % [nModes x nT]
        GFP = sqrt(sum(F.^2, 1));
        [~, frame] = max(GFP);
        Jframe = Res.ImagingKernel * F(:, frame);       % [3nV x 1]
    end

    out = struct('S', S, 'basis', basis, 'dbasis', dbasis, 'Results', Res, ...
                 'c', c, 'F', F, 'Time', rec.Time, 'sfreq', rec.sfreq, ...
                 'iSel', iSel, 'SurfaceFile', cortexFile, 'frame', frame, 'Jframe', Jframe, ...
                 'Comment', rec.Comment, 'Gm', Gm, 'NoiseCov', Cnoise, ...
                 'studyDir', studyDirOut, 'dataName', dataNameOut);
    out.ChannelTypes = chTypes;      % cell -> assign after struct() to avoid a struct array
end

% ----- helpers -----
function f = i_find(d, pat)
    L = dir(fullfile(d, pat));  f = '';
    L = L(~startsWith({L.name}, '._'));                 % drop AppleDouble sidecar files
    if isempty(L), return; end
    [~, newest] = max([L.datenum]);                     % prefer the freshest match (recomputed headmodel)
    f = fullfile(d, L(newest).name);
end
function A = i_sub(A, idx)
    if ~isempty(A) && size(A,1) >= max(idx), A = A(idx, idx); end
end

% Author: Diellor Basha, 2026
