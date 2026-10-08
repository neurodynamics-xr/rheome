function S = scatter(basis, F, frame, opts)
% FILTERS.SCATTER  Graph wavelet scattering: |W| cascaded through the frame, then pooled.
%
%   S = rheome.filters.scatter(basis, F, frame)
%   S = rheome.filters.scatter(basis, F, frame, opts)
%
% The scattering transform (Mallat 2012) is a convolutional network whose filters are
% FIXED to wavelets: convolve, take the modulus, convolve again, average. It has NO
% trainable parameters, which is what makes it usable on a handful of subjects where a
% learned network cannot be fit at all.
%
%   S0(t)       = < F >                              the pooled field
%   S1(m1,t)    = < |g_m1(L) F| >                    one modulus
%   S2(m1,m2,t) = < |g_m2(L) |g_m1(L) F|| >          two
%
% <.> is the mass-weighted spatial average, so every coefficient is a NUMBER per frame,
% not a map: the output is [M] and [M x M] wide rather than [nV]. That collapse is the
% point -- it buys invariance to where on the cortex the structure sat, and it is what
% makes the coefficient count small enough to do statistics on.
%
% ⭐ WHY THE SECOND ORDER IS THE WHOLE REASON TO BUILD THIS. A phase-randomised surrogate
% preserves the power spectrum EXACTLY and sets every higher-order statistic to its
% Gaussian value. S1 is, to leading order, a spectral readout -- pooled band energy --
% so a real recording and its surrogate agree on S1 BY CONSTRUCTION and S1 can never
% discriminate them. S2 measures cross-scale envelope co-modulation, which the power
% spectrum does not determine. Any real-vs-surrogate difference must therefore appear at
% second order or not at all. That makes "S1 matches, S2 matches too" a far stronger null
% than a list of hand-chosen statistics: it is a claim about a whole function class.
%
% ⭐ PASS THE ANALYTIC SIGNAL. On a graph there is no Hilbert transform, so a real mexhat
% coefficient oscillates and |.| RECTIFIES it (doubling the spatial frequency) instead of
% returning an envelope -- the known weak point of graph scattering. Here that is dodged
% for free: F is already complex when it comes from the temporal Hilbert transform
% (rheome.flowpage.complexmaps, rheome.flow.demod), so |.| is a genuine envelope. Feeding a real field
% still runs, but S1 then mixes in the rectification artefact.
%
% ⭐ THE SECOND WAVELET MUST BE COARSER. |g_m1(L) F| is a smooth envelope whose energy
% lies BELOW m1's band, so a finer m2 sees essentially nothing -- Mallat's lambda2 <
% lambda1 rule. Pairs with sigma(m2) <= sigma(m1) are left NaN rather than computed,
% which is also what keeps the second order affordable (roughly half of M^2).
%
% ⭐ S1's MEMBER INDEX IS A CALIBRATED SCALE READOUT -- INSIDE THE BAND. Plant an atom
% built from member m and S1 peaks on m exactly, but only for members between the basis
% floor (frame.SigmaFloor) and roughly 0.6 x the domain radius; outside that the peak
% slides by one, because the basis cannot resolve the fine end and the coarse members
% saturate against the surface's own size. Note it is the POOLED MODULUS that is matched:
% over the same band the L2 energy sum|g_m c|^2 peaks one member FINER, so do not
% substitute a spectral-energy proxy for S1 and expect the same scale. See tScatter.
%
% ⚠ MEMORY IS BOUNDED BY CONSTRUCTION, NOT BY HOPE. The modulus has to be taken at
% VERTICES (|Phi c| ~= Phi |c|), so an [nV x nT] intermediate is unavoidable -- at cortex
% scale that is hundreds of MB per member, and there are M of them. This function
% therefore processes TIME IN BLOCKS sized from opts.MaxBytes and never holds more than
% one member's field at once, so peak memory is set by the budget and NOT by nT. The
% block size is reported in S.BlockSize; the answer is identical for any of them.
%
% INPUTS:
%   basis  struct with .Phi [nV x K], .Lambda [K x 1], .Mass [nV x nV]
%   F      [nV x nT] field, real or complex (complex strongly preferred -- see above)
%   frame  from rheome.filters.frame (needs .g and .Sigma)
%   opts   .Order     1 | 2                        (default 2)
%          .MaxBytes  peak array budget in bytes   (default 5e8 = 500 MB)
%          .Pool      [nV x nR] pooling weights    (default [] = mass-weighted global mean)
%
% OUTPUT (struct S):
%   .S0 [nR x nT]          pooled field (complex if F is)
%   .S1 [M x nR x nT]      first order, real and >= 0
%   .S2 [M x M x nR x nT]  second order, NaN where the coarser-second rule excludes the pair
%   .Sigma [1 x M]         member spatial scale in metres (Inf for a low-pass member)
%   .Pairs [P x 2]         the (m1,m2) pairs actually computed
%   .Order .BlockSize .NumBlocks .PeakBytes .MaxBytes
%
% With the default pooling nR = 1, so squeeze(S.S1) is [M x nT].
%
% See also: rheome.filters.frame, rheome.filters.frame_analysis, rheome.filters.localize, rheome.jointfilterbank
%
% Author: Diellor Basha, 2026 (after Mallat 2012; Bruna & Mallat 2013)

    if nargin < 4 || isempty(opts), opts = struct(); end
    order    = i_opt(opts, 'Order',    2);
    maxBytes = i_opt(opts, 'MaxBytes', 5e8);
    poolW    = i_opt(opts, 'Pool',     []);

    if ~ismember(order, [1 2])
        error('filters:scatter:order', 'Order must be 1 or 2, got %s.', mat2str(order));
    end
    if ~(isscalar(maxBytes) && maxBytes > 0)
        error('filters:scatter:maxBytes', 'MaxBytes must be a positive scalar.');
    end
    if ~isfield(frame, 'Sigma')
        error('filters:scatter:frame', ...
            'frame needs .Sigma to order the cascade. Build it with rheome.filters.frame.');
    end

    Phi = basis.Phi;  Mass = basis.Mass;  lam = double(basis.Lambda(:));
    nV  = size(Phi, 1);  nT = size(F, 2);
    if size(F, 1) ~= nV
        error('filters:scatter:size', ...
            'F has %d rows but the basis has %d.', size(F,1), nV);
    end

    H = rheome.filters.frame_gains(frame, lam);            % [K x M]
    K = size(H, 1);  M = size(H, 2);

    if isempty(poolW)
        w = full(sum(Mass, 2));                     % vertex areas
        w = w / max(sum(w), realmin);               % -> area-weighted mean
    else
        w = poolW;
        if size(w, 1) ~= nV
            error('filters:scatter:pool', ...
                'Pool has %d rows but the basis has %d.', size(w,1), nV);
        end
    end
    nR = size(w, 2);
    wt = w.';                                       % [nR x nV]

    % Low-pass / scaling members carry no scale; they are the COARSEST thing in the bank,
    % so Inf is the ordering-correct reading of the NaN, not a missing value.
    sig = double(frame.Sigma(:)).';
    sig(~isfinite(sig)) = Inf;

    % --- the byte budget -> a block of time columns ------------------------------------
    % Peak concurrent per column, worst case (complex input):
    %   Cb [K x 1] cplx 16   Z [nV x 1] cplx 16   U1 [nV x 1] real 8
    %   C1 [K x 1] real  8   Z2 [nV x 1] real  8
    perCol = 32 * nV + 24 * K;
    blk    = max(1, min(nT, floor(maxBytes / perCol)));
    if blk < 1
        error('filters:scatter:tooLarge', ...
            ['A single time column needs %s (nV = %d, K = %d), over MaxBytes = %s. ' ...
             'Reduce K, use a coarser surface, or raise MaxBytes.'], ...
            i_hum(perCol), nV, K, i_hum(maxBytes));
    end

    S0 = zeros(nR, nT);
    if ~isreal(F), S0 = complex(S0); end
    S1 = zeros(M, nR, nT);
    S2 = [];
    if order >= 2, S2 = nan(M, M, nR, nT); end

    pairs = zeros(0, 2);
    for m1 = 1:M
        for m2 = 1:M
            if sig(m2) > sig(m1), pairs(end+1, :) = [m1 m2]; end %#ok<AGROW>
        end
    end

    nb = 0;
    for b0 = 1:blk:nT
        j  = b0:min(b0 + blk - 1, nT);
        nj = numel(j);
        nb = nb + 1;

        Fb = F(:, j);
        Cb = Phi' * (Mass * Fb);                    % [K x nj]
        S0(:, j) = wt * Fb;
        clear Fb

        for m1 = 1:M
            U1 = abs(Phi * (H(:, m1) .* Cb));       % [nV x nj] real envelope
            S1(m1, :, j) = reshape(wt * U1, 1, nR, nj);

            if order >= 2
                C1 = Phi' * (Mass * U1);            % [K x nj]
                clear U1
                for m2 = 1:M
                    if ~(sig(m2) > sig(m1)), continue; end
                    U2 = abs(Phi * (H(:, m2) .* C1));
                    S2(m1, m2, :, j) = reshape(wt * U2, 1, 1, nR, nj);
                end
                clear C1 U2
            else
                clear U1
            end
        end
        clear Cb
    end

    S = struct('S0', S0, 'S1', S1, 'S2', S2, ...
               'Sigma', sig, 'Pairs', pairs, 'Order', order, ...
               'BlockSize', blk, 'NumBlocks', nb, ...
               'PeakBytes', perCol * blk, 'MaxBytes', maxBytes);
end

function v = i_opt(o, f, d)
    if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function s = i_hum(b)
    if     b >= 1e9, s = sprintf('%.2f GB', b/1e9);
    elseif b >= 1e6, s = sprintf('%.2f MB', b/1e6);
    else,            s = sprintf('%.1f kB', b/1e3);
    end
end

% Author: Diellor Basha, 2026
