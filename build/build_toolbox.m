function file = build_toolbox(outDir)
% BUILD_TOOLBOX  Package the toolbox as rheome-<version>.mltbx (R2023b packageToolbox).
%
%   file = build_toolbox()           % writes build/out/rheome-1.0.1.mltbx
%   file = build_toolbox(outDir)
%
% Reproducible: the files are the GIT-TRACKED files of the packaged folders at HEAD, so untracked
% caches (+data/), results/, .asv files and macOS ._ forks never get in. Refuses a dirty tree, so the
% .mltbx always matches a commit. The version is read from the root Contents.m, the one place it lives.
%
% In the .mltbx: +rheome/ (every package and class), doc/ without its sources and build tools,
% Contents.m, README.md, CITATION.cff and LICENSE. Out: tests/, doc/source/, doc/tools/ and build/ --
% they ship in the git repository (and its Zenodo archive), not in the add-on. doc/ is on the path so
% the Help browser finds doc/info.xml; no search index is built, since it records the build machine's
% paths.
%
% Run from anywhere:  matlab -batch "run('<repo>/build/build_toolbox.m')"  or call it.
%
% See also: check_install, matlab.addons.toolbox.ToolboxOptions, matlab.addons.toolbox.packageToolbox
%
% Author: Diellor Basha, 2026

    root = fileparts(fileparts(mfilename('fullpath')));
    if nargin < 1 || isempty(outDir), outDir = fullfile(root, 'build', 'out'); end

    [rc, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no -- . ":!.claude"', root));
    assert(rc == 0, 'build_toolbox:git', 'git is needed to list the tracked files: %s', dirty);
    assert(isempty(strtrim(dirty)), 'build_toolbox:dirty', ...
        'Commit first; the .mltbx must match a commit. Modified:\n%s', dirty);

    c = fileread(fullfile(root, 'Contents.m'));
    v = regexp(c, '% Version (\d+\.\d+\.\d+)', 'tokens', 'once');
    assert(~isempty(v), 'build_toolbox:version', 'No "%% Version x.y.z" line in Contents.m');
    version = v{1};

    % ⚠ not ls-files -z: system() does not return the NULs intact
    [~, ls] = system(sprintf('git -C "%s" -c core.quotePath=false ls-files', root));
    rel = splitlines(strtrim(ls))';
    top = extractBefore(string(rel) + "/", "/");
    keep = top == "+rheome" ...
         | (top == "doc" & ~startsWith(string(rel), ["doc/source/", "doc/tools/"])) ...
         | ismember(string(rel), ["Contents.m", "README.md", "CITATION.cff", "LICENSE"]);
    files = fullfile(root, rel(keep));

    % ⚠ the identifier is the toolbox's identity across versions: never regenerate it
    opts = matlab.addons.toolbox.ToolboxOptions(root, 'c909aace-d971-4726-8aeb-5c95a2dcb320');
    opts.ToolboxName          = 'Rheome';
    opts.ToolboxVersion       = version;
    opts.Summary              = 'The multiscale geometry of human cortical dynamics: graph wavelets, phase and flow on the cortex.';
    opts.Description          = ['Graph wavelet frames (rheome.graphfilterbank, rheome.jointfilterbank, ' ...
        'rheome.timefilterbank), the joint (lambda, omega) plane, phase geometry and cortical flow from MEG ' ...
        'source estimates, on cortical meshes, connectomes and sensor arrays. Mirrors MATLAB''s ' ...
        'cwtfilterbank API on graphs and meshes. Everything lives in one namespace, rheome.*.'];
    opts.AuthorName           = 'Diellor Basha';
    opts.AuthorCompany        = 'McGill University';
    opts.ToolboxFiles         = files;
    opts.ToolboxMatlabPath    = {root, fullfile(root, 'doc')};   % +rheome resolves from the root; doc/info.xml
    opts.ToolboxGettingStartedGuide = fullfile(root, 'doc', 'GettingStarted.mlx');
    opts.MinimumMatlabRelease = 'R2023b';
    if ~isfolder(outDir), mkdir(outDir); end
    opts.OutputFile           = fullfile(outDir, sprintf('rheome-%s.mltbx', version));

    matlab.addons.toolbox.packageToolbox(opts);
    file = opts.OutputFile;
    [~, head] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
    fprintf('built %s from %s (%d files)\n', file, strtrim(head), numel(files));
end

% Author: Diellor Basha, 2026
