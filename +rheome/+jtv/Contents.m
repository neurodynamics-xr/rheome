% +JTV  Joint time-vertex filterbanks: construction, duals, analysis and synthesis.
%
% A joint filter is a weighting of the joint spectrum C(lambda,omega). This package makes a BANK
% of them invertible, which is what turns joint filtering from a one-way operation into a
% transform: decompose into joint bands, work on each, and put them back together exactly.
%
%   bank       - build a bank from separable and non-separable kernels (the dgw composition)
%   bounds     - frame operator S = sum |g_i|^2, bounds A and B, tightness, coverage
%   dual       - canonical dual gd_i = g_i / S  (perfect reconstruction)
%   analysis   - C -> per-member coefficients
%   synthesis  - per-member coefficients -> C
%
% PERFECT RECONSTRUCTION. With gd the canonical dual, sum_i gd_i * g_i = 1 wherever S > 0, so
% synthesis(analysis(C)) = C exactly. Without a dual the bank is analysis-only and nothing can be
% put back.
%
% Ported from GSPBox's gsp_jtv_* family (Francesco Grassi, 2016) -- gsp_jtv_design_can_dual,
% gsp_jtv_evaluate_can_dual, gsp_jtv_design_dgw -- restated for a [K x Omega] joint spectrum
% rather than a vectorised N*T signal.
%
% ⚠ COVERAGE IS A PRECONDITION. The dual is g_i/S, so it is defined only where S > 0. A bank that
% leaves part of the (lambda,omega) plane uncovered has A = 0 and cannot be inverted there. Always
% read rheome.jtv.bounds before trusting rheome.jtv.dual.
%
% See also: rheome.filters.frame, rheome.flow.joint
%
% Author: Diellor Basha, 2026 (after F. Grassi, GSPBox)

% Author: Diellor Basha, 2026
