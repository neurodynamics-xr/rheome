% IO  Input/output for the atom-design demo.
%
% Reads Brainstorm data from disk into plain MATLAB structs, with no dependency on
% the Brainstorm environment (no bst_get, protocol database, or GUI).
%
% Sub-packages:
%   rheome.io.read.surface   - read a Brainstorm cortical surface .mat
%   rheome.io.read.fssurf    - read a FreeSurfer binary surface (lh.white, lh.sphere, ...)
%   rheome.io.read.atlas     - read every parcellation (the ROI axis) from that surface
%   rheome.io.read.nxrsurface - a subject's cortex and medial wall from an nxr store (the editor's dataset folder)
%   rheome.io.read.channel   - sensor definitions        rheome.io.read.headmodel - leadfield
%   rheome.io.read.recording - sensor time series        rheome.io.read.noisecov  - noise covariance
%   rheome.io.read.inverse   - a kernel-only inverse results file
%   rheome.io.read.fibers    - fiber tracts
%   rheome.io.write.bundle   - a result on a surface as a rheome.bundle@1 folder (what rheome.view opens)
%
% Author: Diellor Basha, 2026
