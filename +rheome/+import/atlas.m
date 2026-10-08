function file = atlas(name)
% IMPORT.ATLAS  Read the surface's parcellations and cache them to +data.
%
%   file = rheome.import.atlas(name)
%
% Reads every atlas out of the dataset's Brainstorm cortical surface (rheome.io.read.atlas) and
% writes them to +data/<name>/atlas.mat, so the ROI axis is available without the external
% volume being mounted. Follows the same one-file-per-concern pattern as rheome.import.bases /
% rheome.import.connectome / rheome.import.sphere.
%
% ⚠ REQUIRES THE ORIGINAL SURFACE FILE. The cached surface.mat keeps geometry and
% operators but NOT the atlases -- rheome.io.read.surface deliberately extracts only the
% 'Structures' hemisphere split. So this reads `cortexFile` from surface.mat and goes back
% to the source. If that volume is unmounted, rheome.import.atlas is the one step that fails.
%
% See also: rheome.io.read.atlas, rheome.load.atlas, rheome.import.surface, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    name  = char(name);
    dsdir = fullfile(rheome.load.root(), name);
    sfile = fullfile(dsdir, 'surface.mat');
    if ~exist(sfile, 'file')
        error('import:atlas:missing', ...
            'No surface cached for ''%s''. Run  rheome.import.surface / rheome.import.dataset  first.', name);
    end

    C = builtin('load', sfile, 'cortexFile');
    if ~isfield(C, 'cortexFile') || isempty(C.cortexFile) || ~exist(C.cortexFile, 'file')
        error('import:atlas:cortexFile', ...
            ['Cached surface for ''%s'' points at a cortex file that is not reachable:\n  %s\n' ...
             'Mount the dataset volume, or re-import the surface from its current path.'], ...
            name, i_str(C));
    end

    fprintf('rheome.import.atlas[%s]: reading parcellations\n', name);
    atlas = rheome.io.read.atlas(C.cortexFile);              %#ok<NASGU>
    cortexFile = C.cortexFile;                        %#ok<NASGU>

    file = fullfile(dsdir, 'atlas.mat');
    builtin('save', file, 'atlas', 'cortexFile', '-v7.3');

    for k = 1:numel(atlas)
        a = atlas(k);
        if a.nScout == 0
            fprintf('  %-22s (empty)\n', a.Name);
        else
            fprintf('  %-22s %3d scouts, %5d/%d vertices (%.1f%%)%s\n', ...
                a.Name, a.nScout, a.nCovered, a.nV, 100*a.Coverage, ...
                i_flag(a.Disjoint));
        end
    end
    fprintf('  -> %s\n', file);
end

function s = i_flag(disjoint)
    if disjoint, s = ''; else, s = '  ** OVERLAPPING **'; end
end

function s = i_str(C)
    if isfield(C, 'cortexFile') && ~isempty(C.cortexFile), s = C.cortexFile; else, s = '(none)'; end
end

% Author: Diellor Basha, 2026
