function R = dirac(Gm, dbasis, NoiseCovMat, opts)
% INVERSE.DIRAC  Whitened minimum-norm inverse in the Dirac eigenmode basis.
%
%   R = rheome.inverse.dirac(Gm, dbasis, NoiseCovMat)
%   R = rheome.inverse.dirac(Gm, dbasis, NoiseCovMat, opts)
%
% Faithful port of Brainstorm's bst_inverse_dirac (which re-implements the whitened
% minimum-norm of bst_inverse_linear_2018 with the source side in the Dirac
% eigenbasis). Five stages:
%   1. Noise-covariance regularization + whitener iW           (reg/diag/none/median/shrink)
%   2. Observability SVD  iW*Gm = UL*S*VL'  + SNR -> Lambda
%   3. Wiener spectral window g(s) = Lambda*s / (Lambda*s^2 + 1)
%   4. Vertex resolution modes  Wres = Phi*VL ; whitened kernel Wres*diag(g)*UL'
%   5. Source normalization (amplitude / dspm2018 / sloreta) + fold in whitener
%
% INPUTS:
%   Gm           [nCh x nModes] mode-forward from rheome.forward.dirac
%   dbasis       Dirac eigenbasis from rheome.eigen.dirac_frame (needed for STAGE 4 reconstruct)
%   NoiseCovMat  struct with .NoiseCov [nCh x nCh] (and .FourthMoment/.nSamples for shrink)
%   opts (all optional):
%     .NoiseMethod   'reg'(default) | 'diag' | 'none' | 'median' | 'shrink'
%     .NoiseReg      ridge fraction for 'reg' (default 0.1)
%     .ChannelTypes  {1 x nCh} type per row (default all 'MEG'); whitened per modality
%     .SnrFixed      assumed SNR (amplitude ratio, default 3 -> SNR^2 = 9)
%     .InverseMeasure 'amplitude'(default) | 'dspm2018' | 'sloreta'
%
% OUTPUT (struct R):
%   .ImagingKernel      [3nV x nCh]   J = ImagingKernel * F  (source 3-vectors from raw data)
%   .ImagingKernelMode  [nModes x nCh] c = ImagingKernelMode * F  (mode coeffs from raw data)
%   .Whitener [nCh x nCh], .UL, .SL, .Lambda, .SNR, .NoiseRankKept, .RankLeadfield
%
% See also: rheome.forward.dirac, rheome.forward.reconstruct
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct(); end
    Def = struct('NoiseMethod','reg','NoiseReg',0.1,'ChannelTypes',[], ...
                 'SnrFixed',3,'InverseMeasure','amplitude');
    opts = i_defaults(opts, Def);

    Gm  = double(Gm);
    nCh = size(Gm, 1);
    if isempty(opts.ChannelTypes), opts.ChannelTypes = repmat({'MEG'}, 1, nCh); end
    C_noise = NoiseCovMat.NoiseCov;
    if size(C_noise,1) ~= nCh
        error('inverse:dirac:size', 'NoiseCov is %dx%d but the mode-forward has %d channels.', ...
            size(C_noise,1), size(C_noise,2), nCh);
    end

    % ===== STAGE 1: noise regularization + whitener (per modality) =====
    Var_noise = diag(C_noise);
    if (norm(C_noise,'fro') - norm(Var_noise,'fro')) < eps(single(norm(Var_noise,'fro')))
        opts.NoiseMethod = 'diag';          % numerically diagonal -> enforce diag (MNE behaviour)
    end
    if strcmpi(opts.NoiseMethod, 'diag'), C_noise = diag(diag(C_noise)); end

    types = unique(opts.ChannelTypes);
    iW_noise = zeros(nCh);
    NoiseRankKept = 0;
    FM = []; nS = [];
    if strcmpi(opts.NoiseMethod,'shrink')
        FM = NoiseCovMat.FourthMoment; nS = NoiseCovMat.nSamples(1);
    end
    for i = 1:numel(types)
        ndx = find(strcmpi(types{i}, opts.ChannelTypes));
        Csub = (C_noise(ndx,ndx) + C_noise(ndx,ndx)')/2;
        if isempty(FM), FMsub = []; else, FMsub = FM(ndx,ndx); end
        iW_noise(ndx,ndx) = whiten(Csub, opts.NoiseMethod, opts.NoiseReg, FMsub, nS);
        NoiseRankKept = NoiseRankKept + sum(abs(diag(iW_noise(ndx,ndx))) > 0);
    end
    GainWhitened = iW_noise * Gm;

    % ===== STAGE 2: observability SVD + SNR -> Lambda =====
    [UL, Ssvd, VL] = svd(GainWhitened, 'econ');
    SL  = diag(Ssvd);
    SL2 = SL.^2;
    RankLeadfield = sum(SL > length(SL)*eps(single(SL(1))));
    SNR    = opts.SnrFixed^2;
    Lambda = SNR / mean(SL2);                          % Hamalainen mean-eigenvalue

    % ===== STAGE 3: Wiener window on the observability axis =====
    gWin       = (Lambda * SL) ./ (Lambda * SL2 + 1);  % [r x 1]
    KernelMode = VL * (gWin .* UL');                   % [nModes x nCh] (whitened data)

    % ===== STAGE 4: vertex resolution modes + current kernel =====
    Wres  = rheome.forward.reconstruct(VL, dbasis);           % [3nV x r]
    nV    = size(Wres, 1) / 3;
    KvtxW = Wres * (gWin .* UL');                       % [3nV x nCh] whitened-data current kernel

    % ===== STAGE 5: source normalization + fold in whitener =====
    switch lower(opts.InverseMeasure)
        case 'amplitude'
            Knorm = KvtxW;
        case 'dspm2018'
            rowN2  = sum(KvtxW.^2, 2);
            Snorm  = sqrt(sum(reshape(rowN2, 3, []), 1))';
            Snorm3 = reshape(repmat(Snorm', 3, 1), [], 1); Snorm3(Snorm3 < eps) = eps;
            Knorm  = KvtxW ./ Snorm3;
        case 'sloreta'
            Lw = iW_noise * rheome.forward.reconstruct(Gm.', dbasis).';   % [nCh x 3nV] whitened effective leadfield
            Knorm = zeros(size(KvtxW));
            for v = 1:nV
                idx = (v-1)*3 + (1:3);
                Rv  = KvtxW(idx,:) * Lw(:,idx);
                [Ur,Sr,Vr] = svd((Rv+Rv')/2); sr = diag(Sr);
                rk = sum(sr > length(sr)*eps(single(sr(1))));
                SIR = Vr(:,1:rk) * diag(1./sqrt(sr(1:rk))) * Ur(:,1:rk)';
                Knorm(idx,:) = SIR * KvtxW(idx,:);
            end
        otherwise
            error('inverse:dirac:measure', 'InverseMeasure must be amplitude/dspm2018/sloreta.');
    end

    R = struct();
    R.ImagingKernel     = Knorm * iW_noise;      % [3nV x nCh] on RAW data
    R.ImagingKernelMode = KernelMode * iW_noise; % [nModes x nCh] on RAW data
    R.Whitener      = iW_noise;
    R.UL = UL;  R.SL = SL;  R.Lambda = Lambda;  R.SNR = SNR;  R.WienerWin = gWin;
    R.NoiseRankKept = NoiseRankKept;
    R.RankLeadfield = RankLeadfield;
    R.InverseMeasure = opts.InverseMeasure;
end

% ---------------------------------------------------------------------------

% ----- shared helpers now live in +inverse/private/whiten.m -----

% Author: Diellor Basha, 2026
