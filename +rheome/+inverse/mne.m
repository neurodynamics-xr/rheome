function R = mne(Gain, NoiseCovMat, opts)
% INVERSE.MNE  Whitened minimum-norm inverse on the RAW unconstrained leadfield.
%
%   R = rheome.inverse.mne(Gain, NoiseCovMat)
%   R = rheome.inverse.mne(Gain, NoiseCovMat, opts)
%
% The standard distributed inverse, with no change of source basis: the estimate is the
% smallest-norm current consistent with the whitened data. This is the PRIMITIVE the flow
% operators are built on -- rheome.inverse.dirac is the alternative, not the baseline.
%
% Same five stages as rheome.inverse.dirac, and deliberately the same code path for four of them:
%   1. Noise-covariance regularization + whitener   (+inverse/private/whiten.m, shared)
%   2. Observability SVD   iW*Gain = UL*S*VL'   + SNR -> Lambda
%   3. Wiener spectral window  g(s) = Lambda*s / (Lambda*s^2 + 1)
%   4. Current kernel  VL*diag(g)*UL'
%   5. Source normalization (amplitude / dspm2018 / sloreta) + fold in the whitener
%
% ⭐ THE ONLY DIFFERENCE FROM rheome.inverse.dirac IS THE SOURCE BASIS, AND THAT IS THE POINT. There,
% stage 2 runs on Gm = rheome.forward.dirac(Gain, dbasis) -- the leadfield projected onto K Dirac
% modes -- and stage 4 reconstructs through the basis. Here stage 2 runs on the leadfield
% itself and stage 4 needs no reconstruction, because VL already lives on the 3nV source
% space. Everything else, including the whitener, is literally the same code. So a difference
% between the two results is attributable to the basis and to nothing else.
%
% ⚠ WHAT THE DIRAC ROUTE COSTS, AND WHY THIS IS THE DEFAULT. Projecting onto K modes
% band-limits the current BEFORE the inverse is solved, and that limit is invisible
% downstream: at the pipeline defaults (400 Dirac modes per hemisphere) the reconstructed
% current stops at ~90 mm, while the LBO analysis axis built on top of it reaches ~28 mm --
% a 3.2x mismatch, so the finest apertures were filtering structure the source estimate could
% not represent. The raw inverse has no such hidden floor: its source space is the full 3nV,
% and the only band limit left is the one the ANALYSIS basis imposes, where it is visible and
% reported (rheome.flow.context prints the frame floor).
%
% ⚠ STILL NOT UNBIASED, JUST UNHIDDEN. Minimum norm has its own prior -- it returns the
% zero-null-space representative, which is why rheome.forward.measurable scores any MNE output 1.0
% by construction (see tMeasurable). Depth bias is untouched here; use rheome.flow.sensitivity as a
% mask and rheome.flow.crosstalk for attribution, exactly as before.
%
% INPUTS:
%   Gain         [nCh x 3nV] UNCONSTRAINED leadfield (rheome.io.read.headmodel), rows = the channels
%                the noise covariance is over
%   NoiseCovMat  struct with .NoiseCov [nCh x nCh] (and .FourthMoment/.nSamples for 'shrink')
%   opts (all optional):
%     .NoiseMethod    'reg'(default) | 'diag' | 'none' | 'median' | 'shrink'
%     .NoiseReg       ridge fraction for 'reg' (default 0.1)
%     .ChannelTypes   {1 x nCh} type per row (default all 'MEG'); whitened per modality
%     .SnrFixed       assumed SNR, amplitude ratio (default 3 -> SNR^2 = 9)
%     .nVert          expected vertex count. ⚠ PASS IT WHEN YOU KNOW IT. Shape alone CANNOT
%                     distinguish an unconstrained [nCh x 3nV] gain from a constrained
%                     [nCh x nV] one whenever nV happens to be divisible by 3 -- the
%                     constrained gain is then silently read as covering nV/3 vertices and
%                     every source lands in the wrong place with no error. rheome.flow.context knows
%                     nV from the surface and always passes it.
%     .InverseMeasure 'dspm2018'(default) | 'amplitude' | 'sloreta'
%                     ⚠ rheome.flow.* wants 'amplitude' -- physical current in A.m. The default
%                     matches rheome.inverse.dirac so the two are swappable without a silent change.
%
% OUTPUT (struct R):
%   .ImagingKernel [3nV x nCh]  J = ImagingKernel * F
%   .Whitener [nCh x nCh] .UL .SL .Lambda .SNR .NoiseRankKept .RankLeadfield .Measure
%
% ⚠ MEMORY. Stage 2 is an economy SVD of an [nCh x 3nV] matrix, so VL is [3nV x nCh] --
% 133 MB at 20484 vertices and 270 channels, in double. The kernel is the same size again.
% That is the whole cost; nothing here is [3nV x 3nV].
%
% See also: rheome.inverse.dirac, rheome.forward.measurable, rheome.flow.context, rheome.flow.sensitivity
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    Def = struct('NoiseMethod','reg','NoiseReg',0.1,'ChannelTypes',[], ...
                 'SnrFixed',3,'InverseMeasure','dspm2018','nVert',[]);
    f = fieldnames(Def);
    for i = 1:numel(f)
        if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = Def.(f{i}); end
    end

    Gain = double(Gain);
    nCh  = size(Gain, 1);
    if mod(size(Gain,2), 3) ~= 0
        error('inverse:mne:constrained', ...
            ['Gain is [%d x %d]; the flow path needs an UNCONSTRAINED leadfield [nCh x 3nV]. ' ...
             'rheome.io.read.headmodel returns one.'], nCh, size(Gain,2));
    end
    nV = size(Gain,2) / 3;
    if ~isempty(opts.nVert) && nV ~= opts.nVert
        error('inverse:mne:constrained', ...
            ['Gain is [%d x %d], i.e. %d vertices unconstrained, but %d were expected. ' ...
             'A CONSTRAINED [nCh x nV] gain with nV divisible by 3 reaches here silently -- ' ...
             'pass the unconstrained leadfield from rheome.io.read.headmodel.'], ...
            nCh, size(Gain,2), nV, opts.nVert);
    end
    if isempty(opts.ChannelTypes), opts.ChannelTypes = repmat({'MEG'}, 1, nCh); end

    C_noise = NoiseCovMat.NoiseCov;
    if ~isequal(size(C_noise), [nCh nCh])
        error('inverse:mne:size', ...
            'NoiseCov is %s but the leadfield has %d channels. Subset them to the same rows.', ...
            mat2str(size(C_noise)), nCh);
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
        ndx  = find(strcmpi(types{i}, opts.ChannelTypes));
        Csub = (C_noise(ndx,ndx) + C_noise(ndx,ndx)')/2;
        if isempty(FM), FMsub = []; else, FMsub = FM(ndx,ndx); end
        iW_noise(ndx,ndx) = whiten(Csub, opts.NoiseMethod, opts.NoiseReg, FMsub, nS);
        NoiseRankKept = NoiseRankKept + sum(abs(diag(iW_noise(ndx,ndx))) > 0);
    end
    GainWhitened = iW_noise * Gain;                    % [nCh x 3nV]

    % ===== STAGE 2: observability SVD + SNR -> Lambda =====
    [UL, Ssvd, VL] = svd(GainWhitened, 'econ');        % VL [3nV x r]
    SL  = diag(Ssvd);
    SL2 = SL.^2;
    RankLeadfield = sum(SL > length(SL)*eps(single(SL(1))));
    SNR    = opts.SnrFixed^2;
    Lambda = SNR / mean(SL2);                          % Hamalainen mean-eigenvalue

    % ===== STAGE 3-4: Wiener window -> current kernel (no reconstruction needed) =====
    gWin  = (Lambda * SL) ./ (Lambda * SL2 + 1);       % [r x 1]
    KvtxW = VL * (gWin .* UL');                        % [3nV x nCh], whitened-data kernel

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
            Lw = GainWhitened;                          % the whitened effective leadfield IS the gain here
            Knorm = zeros(size(KvtxW));
            for v = 1:nV
                idx = (v-1)*3 + (1:3);
                Rv  = KvtxW(idx,:) * Lw(:,idx);
                [Ur,Sr,Vr] = svd((Rv+Rv')/2); sr = diag(Sr);
                rk  = sum(sr > length(sr)*eps(single(sr(1))));
                SIR = Vr(:,1:rk) * diag(1./sqrt(sr(1:rk))) * Ur(:,1:rk)';
                Knorm(idx,:) = SIR * KvtxW(idx,:);
            end
        otherwise
            error('inverse:mne:measure', 'Unknown InverseMeasure "%s".', opts.InverseMeasure);
    end

    R = struct('ImagingKernel', Knorm * iW_noise, 'Whitener', iW_noise, ...
               'UL', UL, 'SL', SL, 'Lambda', Lambda, 'SNR', SNR, ...
               'NoiseRankKept', NoiseRankKept, 'RankLeadfield', RankLeadfield, ...
               'Measure', lower(opts.InverseMeasure), 'nVert', nV);
end

% Author: Diellor Basha, 2026
