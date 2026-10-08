function scal = frame_scalogram(basis, F, frame)
% FILTERS.FRAME_SCALOGRAM  Per-scale energy of a field over time (a spectral-graph scalogram).
%
%   scal = rheome.filters.frame_scalogram(basis, F, frame)
%
% For each frame member m, the wavelet band is W_m = g_m(L) F and its energy at time t is the
% vertex-sum of squared magnitude, E(m,t) = sum_v |W_m(v,t)|^2. Stacking over m gives a
% SCALOGRAM E [M x nT]: which spatial scale carries the field at each instant (the cortical
% analogue of a wavelet scalogram / Morlet time-scale map). For a Dirac (quaternion) basis the
% magnitude is the PHYSICAL current magnitude -- the imaginary 3-vector (x,y,z) norm, not the
% full quaternion. Energy is computed in coefficient space via the Euclidean Gram G = Phi'*Phi
% (matching bst_eigenwavelet's Scalogram exactly) so the [rows x nT x M] field is never
% materialised -- safe for long time extents. Port of bst_eigenwavelet('ScalogramEnergy').
%
% INPUTS:
%   basis  struct with .Phi [rows x K], .Lambda [K x 1], .Mass (B); Dirac basis also has .nVert
%   F      [rows x nT] field (3-vector embedded as a pure quaternion for a Dirac basis)
%   frame  from rheome.filters.frame
%
% OUTPUT (struct scal):
%   .energy  [M x nT]  per-member energy over time (scales x time)
%   .centers [1 x M]   per-member characteristic wavenumber sqrt(lambda) (label the scale axis)
%   .Family            the frame family
%
% See also: rheome.filters.frame, rheome.filters.frame_analysis, rheome.filters.tovec
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('ScalogramEnergy'))

    Phi = basis.Phi;  B = basis.Mass;
    H   = rheome.filters.frame_gains(frame, basis.Lambda);      % [K x M]
    C   = Phi' * (B * F);                                % [K x nT] mode coefficients
    M   = size(H, 2);  nT = size(F, 2);

    isQuat = isfield(basis, 'nVert') && (size(Phi, 1) == 4 * basis.nVert);
    if isQuat
        rowsE = reshape((0:basis.nVert-1)*4 + (2:4)', [], 1);   % imag (x,y,z) rows = current
        PhiE  = Phi(rowsE, :);
    else
        PhiE  = Phi;
    end
    G = PhiE' * PhiE;                                    % [K x K] Euclidean energy Gram

    energy = zeros(M, nT);
    for m = 1:M
        Dm = H(:, m) .* C;                              % [K x nT]
        energy(m, :) = real(sum(conj(Dm) .* (G * Dm), 1));
    end
    scal = struct('energy', energy, 'centers', frame.Centers, 'Family', frame.Family);
end

% Author: Diellor Basha, 2026
