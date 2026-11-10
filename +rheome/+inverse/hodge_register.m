function Bt = hodge_register(B, T, muS, muT)
% INVERSE.HODGE_REGISTER  A subject's fused weights carried onto a template: sensors straight to the template's synthesis weights.
%
%   Bt = rheome.inverse.hodge_register(B, T, muS, muT)    % flow modes carried as flow FIELDS (the stored fsaverage weights')
%   Bt = rheome.inverse.hodge_register(B, T)              % flow modes carried as POTENTIALS
%
% B    [2Ks x C] the subject's fused weights (rheome.inverse.hodge_weights);  T [Kt x Ks] the transfer (rheome.eigen.transfer)
% muS, muT  the two bases' eigenvalues (pkg.mu)
% Bt   [2Kt x C]: [T a; T b] as potentials, or [T a; muT^-1 T muS b] as flow fields.
%
% The weights are linear in the kernel and T is linear in the coefficients, so registering the weights registers every map
% made from them: w_t(t) = Bt x(t). The normal-current weights a are a scalar's coefficients (T a). The flow weights b are
% potentials, and a potential carried as it is keeps the drop between corresponding points, so on a smaller template its flow
% is steeper by the length ratio, a bias that tracks brain size; carried as flow fields of unit strength (|grad psi_k| = mu_k),
% each mode keeps its flow instead.
% Port of nxr-convert eigen/fused.py registered_weights.
%
% See also: rheome.eigen.transfer, rheome.inverse.hodge_weights, rheome.inverse.hodge_synthesis
%
% Author: Diellor Basha, 2026 (port of nxr-convert registered_weights)
arguments, B (:,:) double, T (:,:) double, muS (:,1) double = [], muT (:,1) double = [], end
Ks = size(T, 2);
if size(B, 1) ~= 2 * Ks, error('rheome:hodge_register', 'weights for %d modes, the transfer is from %d', size(B, 1) / 2, Ks); end
if isempty(muS) ~= isempty(muT), error('rheome:hodge_register', 'the flow-field transfer needs both spectra (muS and muT)'); end
if isempty(muS), Bt = [T * B(1:Ks, :); T * B(Ks+1:end, :)]; return, end
Bt = [T * B(1:Ks, :); (T * (muS(1:Ks) .* B(Ks+1:end, :))) ./ muT(1:size(T, 1))];
end
