function report = build_reference(toolboxRoot, varargin)
%BUILD_REFERENCE  Generate a toolbox's function reference and Help-browser table of contents from its help text.
%   REPORT = BUILD_REFERENCE(TOOLBOXROOT) reads the help text of every class (@folder or classdef file) and every
%   package (+folder) function under TOOLBOXROOT, writes one HTML page per item to TOOLBOXROOT/doc/reference/,
%   writes TOOLBOXROOT/doc/helptoc.xml, and builds the documentation search database. Nothing in doc/reference/
%   or helptoc.xml is edited by hand: fix the help text in the code and run this again.
%
%   Name-value options:
%     'Guide'     - N-by-2 cell {target, title}: hand-written pages under doc/ (tutorial, how-to, explanation)
%                   listed before the reference, in order. Default: {'GettingStarted.html','Getting Started'} if
%                   doc/GettingStarted.html exists, else {}.
%     'Packages'  - cellstr of package names to include. Default: every +folder at the root (not nested).
%     'Exclude'   - cellstr of names (package, class or 'pkg.fn') to leave out. Default: {}.
%     'Title'     - toolbox title for the TOC root. Default: <name> in TOOLBOXROOT/info.xml, else the first word
%                   of TOOLBOXROOT/Contents.m, else the folder name.
%     'Commit'    - the commit the reference is generated from (stamped on every page). Default: git HEAD, else ''.
%     'SearchDB'  - build the search database (builddocsearchdb). Default: true.
%
%   REPORT is a struct: pages (cellstr), missingHelp (items with no help text), noH1 (help text whose first line
%   does not start with the item's name in capitals, NAME or PKG.NAME — the H1 line lookfor and Contents.m use), commit.
%   Missing help and missing H1 lines are gaps to request of the code's owners; never patch them in the docs.
%
%   Example:
%     r = build_reference(pwd, 'Guide', {'GettingStarted.html','Getting Started'; 'howto_flow.html','How to ...'});
%     disp(r.missingHelp)
%
%   See also HELP, BUILDDOCSEARCHDB, PUBLISH, RUN_DOC_EXAMPLES.

p = inputParser;
p.addRequired('toolboxRoot', @(x) isfolder(x));
p.addParameter('Guide', [], @(x) isempty(x) || (iscell(x) && size(x,2) == 2));
p.addParameter('Packages', [], @(x) isempty(x) || iscellstr(x) || isstring(x));
p.addParameter('Exclude', {}, @(x) iscellstr(x) || isstring(x));
p.addParameter('Title', '', @(x) ischar(x) || isstring(x));
p.addParameter('Commit', [], @(x) isempty(x) || ischar(x) || isstring(x));
p.addParameter('SearchDB', true, @islogical);
p.parse(toolboxRoot, varargin{:});
o = p.Results;
root = char(o.toolboxRoot);
docDir = fullfile(root, 'doc');
refDir = fullfile(docDir, 'reference');
if ~isfolder(refDir), mkdir(refDir); end
delete(fullfile(refDir, '*.html'));            % generated: start clean so removed items disappear
addpath(root);
cleanup = onCleanup(@() rmpath(root));
exclude = cellstr(o.Exclude);

% ── provenance ──
commit = char(o.Commit);
if isempty(o.Commit)
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', root));
    if st == 0, commit = strtrim(out); else, commit = ''; end
end
title = char(o.Title);
if isempty(title)
    title = contents_title(root);
end

% ── inventory: classes at the root, then packages ──
items = struct('name', {}, 'kind', {}, 'group', {});
d = dir(root);
for k = 1:numel(d)
    n = d(k).name;
    if d(k).isdir && startsWith(n, '@')
        items(end+1) = struct('name', n(2:end), 'kind', 'class', 'group', 'Classes'); %#ok<AGROW>
    elseif ~d(k).isdir && endsWith(n, '.m') && is_classdef(fullfile(root, n))
        items(end+1) = struct('name', n(1:end-2), 'kind', 'class', 'group', 'Classes'); %#ok<AGROW>
    end
end
if isempty(o.Packages)
    pk = d([d.isdir] & startsWith({d.name}, '+'));
    pkgs = cellfun(@(s) s(2:end), {pk.name}, 'UniformOutput', false);
else
    pkgs = cellstr(o.Packages);
end
for k = 1:numel(pkgs)
    if any(strcmp(pkgs{k}, exclude)), continue; end
    mp = meta.package.fromName(pkgs{k});
    if isempty(mp), warning('build_reference:noPackage', 'Package %s not found on the path.', pkgs{k}); continue; end
    fns = setdiff(sort({mp.FunctionList.Name}), {'Contents'}, 'stable');   % Contents.m is the package's help, not a function
    for j = 1:numel(fns)
        items(end+1) = struct('name', [pkgs{k} '.' fns{j}], 'kind', 'function', 'group', pkgs{k}); %#ok<AGROW>
    end
    cls = sort(arrayfun(@(c) c.Name, mp.ClassList, 'UniformOutput', false));
    for j = 1:numel(cls)
        items(end+1) = struct('name', cls{j}, 'kind', 'class', 'group', 'Classes'); %#ok<AGROW>
    end
end
items = items(~ismember({items.name}, exclude));

% ── one page per item, from its help text ──
report = struct('pages', {{}}, 'missingHelp', {{}}, 'noH1', {{}}, 'commit', commit);
for k = 1:numel(items)
    it = items(k);
    txt = help(it.name);
    short = regexprep(it.name, '^.*\.', '');
    if isempty(strtrim(txt)) || ~isempty(regexp(txt, ['^\s*(\S+\.)?' short ' is a (function|class)'], 'once'))   % MATLAB's stand-in for no help
        report.missingHelp{end+1} = it.name;
        txt = '(No help text. This is a gap in the code, reported to its owners.)';
    else
        first = strtrim(strtok(txt, newline));
        if ~endsWith(['.' strtok(first)], ['.' upper(short)])   % H1 line: NAME or PKG.NAME, then a one-line summary
            report.noH1{end+1} = it.name;
        end
    end
    file = [it.name '.html'];
    write_page(fullfile(refDir, file), it, txt, commit, title);
    items(k).file = file; %#ok<AGROW>  (struct field added on first assignment)
    report.pages{end+1} = fullfile('reference', file);
end

% ── reference index and helptoc.xml ──
write_index(fullfile(refDir, 'index.html'), items, commit, title);
guide = o.Guide;
if isempty(guide)
    if isfile(fullfile(docDir, 'GettingStarted.html')), guide = {'GettingStarted.html', 'Getting Started'}; else, guide = cell(0, 2); end
end
write_toc(fullfile(docDir, 'helptoc.xml'), title, guide, items);
if o.SearchDB
    builddocsearchdb(docDir);
end
end

% ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
function t = contents_title(root)
t = '';
fx = fullfile(root, 'info.xml');                 % the name the Help browser shows
if isfile(fx)
    tok = regexp(fileread(fx), '<name>\s*([^<]+?)\s*</name>', 'tokens', 'once');
    if ~isempty(tok), t = tok{1}; return; end
end
f = fullfile(root, 'Contents.m');
if isfile(f)
    l = strtrim(regexprep(fgetl_safe(f), '^%\s*', ''));
    t = strtrim(strtok(l, ' '));
    if ~isempty(t), return; end
end
[~, t] = fileparts(root);
end

function l = fgetl_safe(f)
fid = fopen(f, 'r'); c = onCleanup(@() fclose(fid));
l = fgetl(fid); if ~ischar(l), l = ''; end
end

function tf = is_classdef(f)
txt = fileread(f);
tf = ~isempty(regexp(txt, '^\s*classdef\>', 'once', 'lineanchors'));
end

function s = esc(s)
s = strrep(s, '&', '&amp;'); s = strrep(s, '<', '&lt;'); s = strrep(s, '>', '&gt;'); s = strrep(s, '"', '&quot;');
end

function write_page(file, it, txt, commit, title)
stamp = '';
if ~isempty(commit), stamp = sprintf(' at commit <code>%s</code>', esc(commit)); end
html = sprintf(['<!DOCTYPE html>\n<html><head><meta charset="utf-8"><title>%s</title></head><body>\n' ...
    '<h1>%s</h1>\n<p><em>%s %s</em> &mdash; %s</p>\n<pre>%s</pre>\n' ...
    '<hr><p><small>Generated from the help text%s by build_reference. Do not edit: fix the help text and rebuild.</small></p>\n' ...
    '<p><a href="index.html">%s reference</a></p>\n</body></html>\n'], ...
    esc(it.name), esc(it.name), esc(it.kind), esc(it.group), esc(title), esc(txt), stamp, esc(title));
write_text(file, html);
end

function write_index(file, items, commit, title)
groups = unique({items.group}, 'stable');
body = '';
for g = 1:numel(groups)
    body = [body sprintf('<h2>%s</h2>\n<ul>\n', esc(groups{g}))]; %#ok<AGROW>
    sel = items(strcmp({items.group}, groups{g}));
    for k = 1:numel(sel)
        body = [body sprintf('<li><a href="%s">%s</a> &mdash; %s</li>\n', esc(sel(k).file), esc(sel(k).name), esc(h1(sel(k).name)))]; %#ok<AGROW>
    end
    body = [body sprintf('</ul>\n')]; %#ok<AGROW>
end
stamp = '';
if ~isempty(commit), stamp = sprintf(' at commit <code>%s</code>', esc(commit)); end
write_text(file, sprintf(['<!DOCTYPE html>\n<html><head><meta charset="utf-8"><title>%s reference</title></head><body>\n' ...
    '<h1>%s reference</h1>\n%s<hr><p><small>Generated from the help text%s by build_reference.</small></p>\n</body></html>\n'], ...
    esc(title), esc(title), body, stamp));
end

function s = h1(name)
t = help(name);
s = strtrim(strtok(t, newline));
short = regexprep(name, '^.*\.', '');
s = strtrim(regexprep(s, ['^' regexptranslate('escape', upper(short)) '\s*'], ''));
end

function write_toc(file, title, guide, items)
x = sprintf('<?xml version="1.0" encoding="utf-8"?>\n<toc version="2.0">\n<tocitem target="%s">%s\n', ...
    esc(guide_first(guide)), esc(title));
for k = 1:size(guide, 1)
    x = [x sprintf('  <tocitem target="%s">%s</tocitem>\n', esc(guide{k,1}), esc(guide{k,2}))]; %#ok<AGROW>
end
x = [x sprintf('  <tocitem target="reference/index.html">Functions and Classes\n')];
groups = unique({items.group}, 'stable');
for g = 1:numel(groups)
    x = [x sprintf('    <tocitem target="reference/index.html">%s\n', esc(groups{g}))]; %#ok<AGROW>
    sel = items(strcmp({items.group}, groups{g}));
    for k = 1:numel(sel)
        x = [x sprintf('      <tocitem target="reference/%s">%s</tocitem>\n', esc(sel(k).file), esc(sel(k).name))]; %#ok<AGROW>
    end
    x = [x sprintf('    </tocitem>\n')]; %#ok<AGROW>
end
x = [x sprintf('  </tocitem>\n</tocitem>\n</toc>\n')];
write_text(file, x);
end

function t = guide_first(guide)
if isempty(guide), t = 'reference/index.html'; else, t = guide{1,1}; end
end

function write_text(file, s)
fid = fopen(file, 'w', 'n', 'UTF-8'); c = onCleanup(@() fclose(fid));
fwrite(fid, s, 'char');
end
