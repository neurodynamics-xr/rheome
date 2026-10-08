% +GRAPHTRANSFORM  Field <-> coefficient transforms for graphfilterbank.
%
% graphfilterbank owns the spectral range and the filter shapes. It does NOT own the
% transform -- exactly as cwtfilterbank does not own fft. This package supplies it.
%
% Everything operator-specific lives behind four handles, which is why Laplace-Beltrami,
% the connection Laplacian and Dirac need no branching anywhere in graphfilterbank:
%   .forward(F)  field -> coefficients      .inverse(C)  coefficients -> field
%   .norm(F)     per-vertex magnitude       .lambda      the exact spectrum (optional)
%
%   rheome.graphtransform.eigen     - exact eigenbasis (Phi, B, Lambda)
%   rheome.graphtransform.chebyshev - polynomial in L; no eigenvectors at all
%   rheome.graphtransform.validate  - shape check
%
% Author: Diellor Basha, 2026
