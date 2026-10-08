function report = build_docs(varargin)
%BUILD_DOCS  Build the toolbox documentation in doc/: the pages, the function reference and helptoc.xml.
%   REPORT = BUILD_DOCS() renders every hand-written page and regenerates the reference:
%     doc/source/GettingStarted.m  -> doc/GettingStarted.mlx and doc/GettingStarted.html (run)
%     doc/examples/howto_*.m       -> doc/howto_*.html (run)
%     doc/source/about_*.m         -> doc/about_*.html (not run: explanation, no code)
%     the help text of every class and package function -> doc/reference/*.html, doc/helptoc.xml
%   Pages are written as live scripts and exported, so figures are embedded in the HTML. The token
%   @COMMIT@ in a page source is replaced by the release tag, v<version> from the root Contents.m.
%   REPORT is build_reference's report (missingHelp, noH1, commit) plus .pages, the HTML written.
%
%   Name-value options:
%     'Pages'      - render the pages (default true)
%     'Reference'  - regenerate the reference and helptoc.xml (default true)
%
%   Run it from a shell, with a display (export needs the Live Editor) and with neutral cache, output
%   and temporary folders, so that no local path is printed into a page:
%     TMPDIR=/tmp/cfdocs/tmp matlab -batch "setenv('RHEOME_DATA','/tmp/cfdocs/data'); \
%        setenv('RHEOME_OUT','/tmp/cfdocs/results'); addpath('doc/tools'); r = build_docs()"
%   Then check every example: t = run_doc_examples(pwd); assert(all(t.Passed)).
%   The search database is not built: builddocsearchdb records the build machine's path.
%
%   See also BUILD_REFERENCE, RUN_DOC_EXAMPLES, EXPORT.

% Author: Diellor Basha, 2026
p = inputParser;
p.addParameter('Pages', true, @islogical);
p.addParameter('Reference', true, @islogical);
p.parse(varargin{:});
o = p.Results;
docDir = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(docDir);
addpath(root, fullfile(docDir, 'examples'));
% stamped with the release (the v<version> tag), not a hash: the published history is rewritten once
v = regexp(fileread(fullfile(root, 'Contents.m')), '% Version (\d+\.\d+\.\d+)', 'tokens', 'once');
commit = ['v' v{1}];

% the guide, in reading order: tutorial, how-to guides, explanations
guide = {
    'GettingStarted.html',          'Getting Started'
    'howto_ingest_brainstorm.html', 'How to import a Brainstorm study'
    'howto_flow_maps.html',         'How to compute flow maps'
    'howto_scale_run.html',         'How to run the per-participant measures'
    'howto_figures.html',           'How to make figures'
    'about_sensor_to_cortex.html',  'From sensors to cortex'
    'about_helmholtz_hodge.html',   'The Helmholtz-Hodge split'
    'about_optical_flow.html',      'Surface optical flow'
    'about_scale.html',             'Scale: local to global, fast to slow'
    'about_gauge.html',             'The gauge'
    'about_resolution_floor.html',  'The resolution floor'
    };

pages = {};
if o.Pages
    tmp = tempname;  mkdir(tmp);  cl = onCleanup(@() rmdir(tmp, 's'));
    src = [dir(fullfile(docDir, 'source', '*.m')); dir(fullfile(docDir, 'examples', 'howto_*.m'))];
    src = src(~startsWith({src.name}, '._'));            % macOS AppleDouble files on external volumes
    for k = 1:numel(src)
        [~, stem] = fileparts(src(k).name);
        txt = strrep(fileread(fullfile(src(k).folder, src(k).name)), '@COMMIT@', commit);
        m = fullfile(tmp, [stem '.m']);  mlx = fullfile(tmp, [stem '.mlx']);
        fid = fopen(m, 'w', 'n', 'UTF-8');  fwrite(fid, txt, 'char');  fclose(fid);
        matlab.internal.liveeditor.openAndSave(m, mlx);      % internal API, works on R2023b
        html = fullfile(docDir, [stem '.html']);
        export(mlx, html, Run=~startsWith(stem, 'about_'));
        close all
        if strcmp(stem, 'GettingStarted'), copyfile(mlx, fullfile(docDir, 'GettingStarted.mlx')); end
        pages{end+1} = html; %#ok<AGROW>
        fprintf('build_docs: %s\n', html);
    end
end

report = struct('pages', {pages});
if o.Reference
    % one namespace: the classes of +rheome, then each of its packages
    pk = dir(fullfile(root, '+rheome', '+*'));
    pkgs = [{'rheome'}, strcat('rheome.', erase({pk.name}, '+'))];
    r = build_reference(root, 'Guide', guide, 'Commit', commit, 'SearchDB', false, 'Title', 'Rheome', ...
        'Packages', pkgs);
    report = r;  report.pages = [pages r.pages];
end
end
% Author: Diellor Basha, 2026
