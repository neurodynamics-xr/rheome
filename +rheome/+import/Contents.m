% IMPORT  Brainstorm -> +data cache (write once). Each rheome.import.* has a rheome.load.* reader.
%
% Functions:
%   rheome.import.atlas        - Read the surface's parcellations and cache them to +data.
%   rheome.import.bases        - Precompute + cache PER-HEMISPHERE operators & eigenbases for both hemispheres.
%   rheome.import.connectome   - Build + cache the whole-brain structural-connectome operators & eigenbases.
%   rheome.import.cortex       - Equip a subject's FOLDED cortical hemisphere with every operator + eigenmode.
%   rheome.import.dataset      - Import a full dataset (surface + study + Dirac eigenbasis) in one call.
%   rheome.import.dirac        - Compute + cache the relative-Dirac eigenbasis for a dataset.
%   rheome.import.inverse      - Find a Brainstorm kernel-only inverse in a study + cache it to +data.
%   rheome.import.rawrecording - Native-rate BST-BIN -> the pageable v7.3 store.
%   rheome.import.recording    - Surface a cached recording as a v7.3, chunk-readable store.
%   rheome.import.sphere       - Equip a subject's FreeSurfer registration sphere with every operator + eigenmode.
%   rheome.import.study        - Read the Brainstorm study structures + cache them to +data.
%   rheome.import.surface      - Read a Brainstorm cortex + compute LBO operators, cache to +data.
%
% Author: Diellor Basha, 2026
