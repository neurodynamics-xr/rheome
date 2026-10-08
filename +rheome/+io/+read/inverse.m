function inv = inverse(resultsFile)
% IO.READ.INVERSE  Read a Brainstorm kernel-only inverse result (results_*_KERNEL_*.mat).
%
%   inv = rheome.io.read.inverse(resultsFile)
%
% Returns the imaging kernel + the metadata needed to reproduce a reconstruction as
%   J = inv.ImagingKernel * Data(inv.GoodChannel, :)
% Kernel-only results only (the file must carry ImagingKernel). Mirrors Brainstorm's
% results struct written by bst_inverse_linear_2018.
%
% OUTPUT (struct inv):
%   .ImagingKernel [nSrc x nGood]      the inverse operator (nSrc = nComponents * nVert)
%   .GoodChannel   [nGood x 1]         data rows the kernel multiplies (indices into the channel file)
%   .nComponents   1 (constrained, normal-oriented) | 3 (unconstrained, 3-vector per vertex)
%   .Function      'mn' | 'dspm2018' | 'sloreta' | ...   (sets the physical units -- see utils.measure_units)
%   .Comment .SurfaceFile .nVert .File
%
% See also: rheome.import.inverse, rheome.load.inverse
%
% Author: Diellor Basha, 2026

    if ~exist(resultsFile, 'file')
        error('io:read:inverse:missing', 'No such results file: %s', resultsFile);
    end
    R = builtin('load', resultsFile);
    if ~isfield(R, 'ImagingKernel') || isempty(R.ImagingKernel)
        error('io:read:inverse:notKernel', ...
            '%s has no ImagingKernel -- only kernel-only results are supported (not full source maps).', resultsFile);
    end
    inv = struct();
    inv.ImagingKernel = R.ImagingKernel;
    inv.GoodChannel   = R.GoodChannel(:);
    inv.nComponents   = double(R.nComponents);
    inv.Function      = char(R.Function);
    inv.Comment       = char(R.Comment);
    inv.SurfaceFile   = i_opt(R, 'SurfaceFile', '');
    inv.nVert         = size(R.ImagingKernel, 1) / max(inv.nComponents, 1);
    inv.File          = resultsFile;
end

function v = i_opt(s, f, d), if isfield(s, f), v = s.(f); else, v = d; end, end

% Author: Diellor Basha, 2026
