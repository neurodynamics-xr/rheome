% LOAD  +data cache -> analyses (read many), and where analysis outputs go (rheome.load.outpath).
%
% Functions:
%   rheome.load.atlas      - Load the cached parcellations (the ROI axis) from +data.
%   rheome.load.bases      - Load the per-hemisphere operators & eigenbases (all families, both hemispheres).
%   rheome.load.connectome - Load the cached whole-brain connectome operators & eigenbases.
%   rheome.load.cortex     - Load the bundled folded-cortex sandbox (operators + eigenmodes).
%   rheome.load.dataset    - Load everything cached for a dataset in one struct.
%   rheome.load.dirac      - Load a cached relative-Dirac eigenbasis (and Laplace–Beltrami basis) from +data.
%   rheome.load.has        - True if 'name' is a cached dataset (a folder under +data), not a file path.
%   rheome.load.inverse    - Load a cached Brainstorm inverse kernel from +data ([] if none).
%   rheome.load.list       - Print the cached datasets in +data and what each holds.
%   rheome.load.outpath    - Where an analysis output goes: <outroot>/[<subject>/]<analysis>/<file>, folder created.
%   rheome.load.outroot    - Absolute path to the ANALYSIS OUTPUT workspace: results, figures, reports.
%   rheome.load.root       - Absolute path to the +data cache folder (shared by rheome.import.* and rheome.load.*).
%   rheome.load.sphere     - Load the bundled registration-sphere sandbox (operators + eigenmodes).
%   rheome.load.study      - Load the cached Brainstorm study structures from +data.
%   rheome.load.surface    - Load a cached cortex surface (and its LBO operators) from +data.
%
% Author: Diellor Basha, 2026
