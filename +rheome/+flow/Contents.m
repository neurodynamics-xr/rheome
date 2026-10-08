% +FLOW  Fused flow imaging kernels (sensors -> cortical flow maps).
%
% Precompute-once, matvec-per-frame operators that map MEG sensor data directly to
% cortical-flow quantities, extending the Dirac ImagingKernelMode fusion through the
% differential (Helmholtz-Hodge) operators. Build the shared context once, then each
% dedicated builder returns its kernel.
%
%   context     - one-time precompute bundle (Dirac source map + currentKernel + operators + LBO)
%   field       - root kernel: current J            [3V x C]   J(:,t)   = vertexOperator * F(:,t)
%   divergence  - sources / sinks                   [V x C]    div(:,t) = vertexOperator * F(:,t)
%   curl        - scalar vorticity                  [V x C]    vort(:,t)= vertexOperator * F(:,t)
%   potential   - Helmholtz Phi (coeff + vertex)    [Ks x C]/[V x C]  sources/sinks potential
%   stream      - Helmholtz Psi (coeff + vertex)    [Ks x C]/[V x C]  stream function (vortices)
%   energy      - total |J|^2 Gram                  [C x C]    E(t)  = F' Q F
%   enstrophy   - total vorticity^2 Gram            [C x C]    Om(t) = F' Q F
%   helicity    - total J.(grad x J) Gram           [C x C]    H(t)  = F' Q F
%   build       - assemble all + coeff-vs-vertex report
%   weak        - internal: Galerkin Bdiv/Brot operators (used by potential/stream)
%   curlvec     - internal: ambient curl-vector kernel (used by helicity)
%
% ⭐⭐ THERE ARE TWO FLOW PIPELINES AND THEY ANSWER DIFFERENT QUESTIONS. Everything above is
% pipeline 1: curl and divergence of the CURRENT, instantaneous, linear, exactly fused, one
% frame in and one map out. Pipeline 2 is kinematic -- how an amplitude PATTERN moves:
%
%   activation  - current -> a scalar activation map [nV x nT] (analytic envelope, or |J|)
%   cortexfeatures - the alpha envelope as mergeable sums per cortical node x store tile
%                 (area, on-area, energy, energy per cortical spatial octave) -> <store>__cortex.mat
%   labelalpha  - alpha class, episodes and trajectories as rheome.select.measure notes on cortical nodes
%   alphastates - query them: one class's states in one hemisphere, joined to their episodes on (level, k)
%   apparent    - that map -> Horn-Schunck velocity v, then div v, curl v, |v|
%                 (reuses rheome.dynamics.opticalflow_scalar and rheome.dynamics.flow_readout)
%
% A rotating current can sit inside a stationary envelope (curl J large, curl v zero) and a
% travelling bump can be irrotational everywhere (curl J zero, div v large). MEASURED on one
% resting alpha tile: the spatial correlation of |curl v| with |curl J| is -0.11, i.e. they are
% unrelated, which is the expected result and not a discrepancy. The vortex-over-the-alpha-
% phase question is pipeline 2. Script: alpha_apparent_flow_omega.m.
%
% ⚠ PIPELINE 2 IS NONLINEAR AND NEEDS TWO FRAMES; pipeline 1 is linear and needs one. Nothing
% about pipeline 2 can be fused into a kernel, which is why it costs 53 s for a 600-frame tile
% on one hemisphere where pipeline 1 costs a matrix-vector product.
%
% Linear kernels are exact fusions of reconstruct + rheome.differential.* (equivalence ~1e-12);
% Grams reproduce the direct sums to ~1e-14. Everything is band-limited by the MNE inverse's
% ~C-dimensional resolution -- fusion is efficiency, not new information.
%
% Author: Diellor Basha, 2026
