function v = vendor_viewer(src)
% VENDOR_VIEWER  Ship the Cortical Flow viewer in Rheome: copy its build into <rheome>/viewer, stamped with its commit.
%
%   v = vendor_viewer()      % from $NXR_VIEWER_DIST, else ~/workspace/research/code/nxr-cortical-flow/editor/dist-viewer
%   v = vendor_viewer(src)
%
% The viewer is Rheome's renderer of the cortex (rheome.view, rheome.snapshot; Diellor, 2026-10-14), so it is versioned
% with Rheome: viewer/ is tracked, packaged by build_toolbox and installed by rheome.view into prefdir/rheome-viewer
% when its VERSION differs. Build the viewer first (npm --prefix editor run build:viewer in nxr-cortical-flow), then run
% this and commit viewer/. The demo bundles (bundles/) are not shipped. Refuses a source repo with uncommitted changes
% under editor/, so VERSION names exactly what was built.
%
% Author: Diellor Basha, 2026
    if nargin < 1 || isempty(src), src = getenv('NXR_VIEWER_DIST'); end
    if isempty(src), src = fullfile(getenv('HOME'), 'workspace', 'research', 'code', 'nxr-cortical-flow', 'editor', 'dist-viewer'); end
    assert(isfile(fullfile(src, 'viewer.html')), 'vendor_viewer:src', 'No viewer build at %s.', src);
    repo = fileparts(fileparts(src));
    [rc, sha] = system(sprintf('git -C "%s" rev-parse --short HEAD', repo));
    assert(rc == 0, 'vendor_viewer:git', 'The viewer build is not in a git repo: %s', repo);
    [~, dirty] = system(sprintf('git -C "%s" status --porcelain -- editor/src editor/package.json', repo));
    assert(isempty(strtrim(dirty)), 'vendor_viewer:dirty', 'Uncommitted viewer sources in %s:\n%s', repo, dirty);

    dst = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'viewer');
    if isfolder(dst), rmdir(dst, 's'); end
    mkdir(dst);
    for f = setdiff({dir(src).name}, {'.', '..', 'bundles', '.DS_Store'})
        copyfile(fullfile(src, f{1}), fullfile(dst, f{1}));
    end
    v = sprintf('nxr-cortical-flow %s', strtrim(sha));
    fid = fopen(fullfile(dst, 'VERSION'), 'w');  fprintf(fid, '%s\n', v);  fclose(fid);
    fprintf('viewer/ <- %s (%s)\n', src, v);
end
