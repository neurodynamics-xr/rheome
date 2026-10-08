function results = run_doc_examples(toolboxRoot, varargin)
%RUN_DOC_EXAMPLES  Run every documentation example of a toolbox and report which pass.
%   RESULTS = RUN_DOC_EXAMPLES(TOOLBOXROOT) runs, each in a fresh workspace, every example script under
%   TOOLBOXROOT/doc/examples/ (*.m), every live script under TOOLBOXROOT/doc/ (*.mlx, GettingStarted.mlx included),
%   and the 'Example:' block of the help text of every item listed with 'HelpOf'. It returns a table with one row
%   per example: Source, Passed, Seconds, Message. Documentation is not done until every row passes.
%
%   Name-value options:
%     'HelpOf'   - cellstr of function or class names whose help-text 'Example:' block is run. Default: {}.
%     'Exclude'  - cellstr of file names to skip (e.g. examples that need private data; they must then say so
%                  and the user-facing page must not depend on them). Default: {}.
%     'Export'   - also export each passing .mlx to HTML next to it (export, R2022a+). Default: false.
%
%   Example:
%     r = run_doc_examples(pwd, 'HelpOf', {'flow.page'});
%     assert(all(r.Passed), 'documentation examples failed')
%
%   See also BUILD_REFERENCE, PUBLISH, EXPORT.

p = inputParser;
p.addRequired('toolboxRoot', @isfolder);
p.addParameter('HelpOf', {}, @(x) iscellstr(x) || isstring(x));
p.addParameter('Exclude', {}, @(x) iscellstr(x) || isstring(x));
p.addParameter('Export', false, @islogical);
p.parse(toolboxRoot, varargin{:});
o = p.Results;
root = char(o.toolboxRoot);
addpath(root); c1 = onCleanup(@() rmpath(root));
docDir = fullfile(root, 'doc');

files = [dir(fullfile(docDir, 'examples', '*.m')); dir(fullfile(docDir, '*.mlx')); dir(fullfile(docDir, '**', '*.mlx'))];
files = files(~ismember({files.name}, cellstr(o.Exclude)));
[~, ia] = unique(fullfile({files.folder}, {files.name}), 'stable');
files = files(ia);

Source = strings(0,1); Passed = false(0,1); Seconds = zeros(0,1); Message = strings(0,1);
for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    [ok, sec, msg] = run_isolated(@() run_script(f));
    Source(end+1,1) = string(erase(f, [root filesep])); Passed(end+1,1) = ok; Seconds(end+1,1) = sec; Message(end+1,1) = msg; %#ok<AGROW>
    if ok && o.Export && endsWith(f, '.mlx')
        export(f, replace(f, '.mlx', '.html'));
    end
end
names = cellstr(o.HelpOf);
for k = 1:numel(names)
    code = help_example(names{k});
    if isempty(code)
        ok = false; sec = 0; msg = "no 'Example:' block in the help text (a gap for the code's owners)";
    else
        [ok, sec, msg] = run_isolated(@() evalin_fresh(code));
    end
    Source(end+1,1) = "help " + names{k}; Passed(end+1,1) = ok; Seconds(end+1,1) = sec; Message(end+1,1) = msg; %#ok<AGROW>
end
results = table(Source, Passed, Seconds, Message);
end

function [ok, sec, msg] = run_isolated(fn)
t = tic; ok = true; msg = "";
fig0 = findall(groot, 'Type', 'figure');
try
    evalc('fn()');
catch err
    ok = false; msg = string(err.message);
end
sec = toc(t);
close(setdiff(findall(groot, 'Type', 'figure'), fig0));
end

function run_script(f)
run(f);       % a script runs in this function's (dynamic) workspace: fresh for each example
end

function evalin_fresh(code)
eval(code);   % runs in this function's own workspace: nothing leaks between examples
end

function code = help_example(name)
% The lines after an 'Example:' or 'Examples:' heading, up to the first blank line or 'See also'.
txt = splitlines(string(help(name)));
i = find(~cellfun(@isempty, regexp(cellstr(txt), '^\s*Examples?:?\s*$', 'once')), 1);
code = '';
if isempty(i), return; end
body = strings(0,1);
for k = i+1:numel(txt)
    l = txt(k);
    if strtrim(l) == "" || startsWith(strtrim(l), "See also"), break; end
    body(end+1) = strtrim(l); %#ok<AGROW>
end
code = char(strjoin(body, newline));
end
