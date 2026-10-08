function W = frame_analysis(basis, F, frame)
% FILTERS.FRAME_ANALYSIS  Wavelet-frame analysis: a field -> its multi-scale coefficients.
%
%   W = rheome.filters.frame_analysis(basis, F, frame)
%
% Projects F onto the eigenbasis once, then filters by each frame member:
%   C = Phi'*(Mass*F);  W(:,:,m) = Phi*(g_m(lambda) .* C)   (= g_m(L) F, member m)
% so W(:,:,m) is F seen through wavelet band m -- the spectral-graph wavelet transform.
% Reconstruct with rheome.filters.frame_synthesis. Port of bst_eigenwavelet('Analysis') for a
% single (per-hemisphere) B-orthonormal basis. Works for scalar (LBO), complex (connection)
% and quaternion (Dirac) bases -- F must already be in that basis's row layout (e.g. a
% [4*nV x nT] quaternion field for Dirac; embed a 3-vector with w = 0 first).
%
% INPUTS:
%   basis  struct with .Phi [rows x K], .Lambda [K x 1], .Mass [rows x rows] (B, mass matrix)
%   F      [rows x nT] field on this basis
%   frame  from rheome.filters.frame
%
% OUTPUT:
%   W      [rows x nT x M] wavelet coefficients, one page per frame member (M = numel(frame.g))
%
% See also: rheome.filters.frame, rheome.filters.frame_synthesis, rheome.filters.apply
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('Analysis'))

    Phi = basis.Phi;  B = basis.Mass;
    H   = rheome.filters.frame_gains(frame, basis.Lambda);      % [K x M]
    C   = Phi' * (B * F);                                % [K x nT]
    M   = size(H, 2);  nT = size(F, 2);
    W   = zeros(size(F, 1), nT, M);
    if ~isreal(F) || ~isreal(Phi), W = complex(W); end
    for m = 1:M
        W(:, :, m) = Phi * (H(:, m) .* C);
    end
end

% Author: Diellor Basha, 2026
