% INVERSE  Source-mapping inverse solutions, and what they can actually resolve.
%
% Whitened minimum-norm estimation. Pure-MATLAB port of Brainstorm's bst_inverse_dirac and
% bst_inverse_linear_2018, with the noise covariance entering as a METRIC (the whitener,
% +inverse/private/whiten.m) rather than as a threshold, and regularisation as a Wiener
% window on the whitened leadfield's singular values.
%
% Functions:
%   rheome.inverse.mne        - whitened MNE on the RAW unconstrained leadfield; the primitive the
%                        flow operators are built on (amplitude / dSPM / sLORETA)
%   rheome.inverse.dirac      - the same five stages with the source side in the Dirac eigenbasis;
%                        a legitimate alternative, not the baseline (rheome.flow.context)
%   rheome.inverse.resolution - the resolution matrix R = K*Gain reduced to MILLIMETRES: PSF spread
%                        by geodesic distance, and the gain per cortical wavelength
%
% ⭐ WHERE THE NOISE FLOOR ACTUALLY LIVES. Three different things get called one:
%   the WHITENER      C^(-1/2), which makes "distance from the data" mean standard deviations
%                     of real noise -- a matrix because MEG noise is spatially correlated
%   the WIENER WINDOW g(s) = s/(s^2 + lambda) with lambda = mean(s^2)/SnrFixed^2, which is a
%                     soft cutoff on the leadfield spectrum: this is the floor MNE computes,
%                     and it is a RANK, not a length
%   the SOURCE NORM   dSPM's k*C*k' per vertex, the variance the kernel yields from noise
%                     alone, which is a floor per LOCATION (and is rheome.flow.sensitivity)
% None of them is a spatial resolution in millimetres. rheome.inverse.resolution is the conversion,
% and it reports that trace(R) is exactly the sum of the Wiener gains -- so the rank the
% regularisation sets and the resolution operator are one object.
%
% Scripts: inverse_resolution_omega.m sweeps SnrFixed against the empty-room band SNRs of
% noise_floor_omega.m and prints both families of length (numbers in
% docs/2026-09-22-ingest-notes.md).
%
% Author: Diellor Basha, 2026
