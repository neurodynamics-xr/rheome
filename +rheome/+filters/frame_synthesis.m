function Frec = frame_synthesis(basis, W, frame, mode)
% FILTERS.FRAME_SYNTHESIS  Wavelet-frame synthesis: multi-scale coefficients -> field.
%
%   Frec = rheome.filters.frame_synthesis(basis, W, frame)          % 'tight' (default)
%   Frec = rheome.filters.frame_synthesis(basis, W, frame, 'dual')  % exact inverse for ANY bank
%
% Filters each wavelet band by a per-member gain and sums:
%   Frec = sum_m Phi*(w_m(lambda) .* (Phi'*(Mass*W(:,:,m))))
%
%   mode = 'tight' (default): w_m = g_m -- the frame operator. Composed with
%       rheome.filters.frame_analysis it applies S(L) = sum_m g_m(lambda)^2, so
%       frame_synthesis(frame_analysis(F)) = S(L) F == A*F when the frame is TIGHT.
%       For 'itersine' A ~ 1 -> EXACT reconstruction (machine precision). This is the exact
%       port of bst_eigenwavelet('Synthesis').
%   mode = 'dual': w_m = g_m / S(L) -- the canonical dual frame. Then
%       frame_synthesis(frame_analysis(F)) = F EXACTLY for ANY bank (itersine|mexhat|heat),
%       wherever the frame covers the spectrum (S > 0). Use this to treat a non-tight
%       scale-space (mexhat/heat) as an invertible wavelet transform.
%
% Single (per-hemisphere) B-orthonormal basis; scalar / complex / quaternion layouts.
%
% INPUTS:
%   basis  struct with .Phi, .Lambda, .Mass  (same basis used for analysis)
%   W      [rows x nT x M] wavelet coefficients from rheome.filters.frame_analysis
%   frame  from rheome.filters.frame;  mode  'tight' (default) | 'dual'
%
% OUTPUT:
%   Frec   [rows x nT] reconstructed field
%
% See also: rheome.filters.frame, rheome.filters.frame_analysis, rheome.filters.frame_bounds
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('Synthesis'))

    if nargin < 4 || isempty(mode), mode = 'tight'; end
    Phi = basis.Phi;  B = basis.Mass;
    H   = rheome.filters.frame_gains(frame, basis.Lambda);      % [K x M]
    switch lower(mode)
        case 'tight'
            Wg = H;
        case 'dual'
            S  = sum(H.^2, 2);                            % frame operator diagonal
            Wg = H ./ max(S, eps);                        % canonical dual gains
        otherwise
            error('filters:frame_synthesis:mode', 'mode must be ''tight'' or ''dual''.');
    end
    M   = size(H, 2);  nT = size(W, 2);
    Frec = zeros(size(W, 1), nT);
    if ~isreal(W) || ~isreal(Phi), Frec = complex(Frec); end
    for m = 1:M
        Frec = Frec + Phi * (Wg(:, m) .* (Phi' * (B * W(:, :, m))));
    end
end

% Author: Diellor Basha, 2026
