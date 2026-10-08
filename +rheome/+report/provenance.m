function P = provenance(opts)
% REPORT.PROVENANCE  Everything needed to reproduce a result, collected rather than typed.
%
%   P = rheome.report.provenance()
%   P = rheome.report.provenance(Script="resolution_scales_omega.m", Data=["...bases.mat"])
%   P = rheome.report.provenance(..., Hash="sha256")     % read every data file and digest it
%
% ⭐ PROVENANCE THAT IS TYPED BY HAND IS PROVENANCE THAT IS WRONG. The commit a figure was
% made from, whether the tree was dirty at the time, and which cached file supplied the data
% are all knowable at run time, so none of them belongs in prose. This returns them as a value
% that rheome.report.article renders into a chapter.
%
% ⚠ A DIRTY TREE IS RECORDED, NOT HIDDEN. .clean is false whenever `git status --porcelain`
% has any output, and .dirtyFiles lists what. rheome.report.article prints a banner in the document
% when this happens, because a result built from uncommitted code cannot be reproduced from
% the SHA alone and a reader has no other way to find that out.
%
% ⚠ HASHING IS OFF BY DEFAULT and that is deliberate on this project: the cached bases are
% 297 MB and digesting them costs seconds per report. Hash="fingerprint" (the default) records
% size and modification time, which detects a changed file but does not prove two files are
% identical. Hash="sha256" proves it. Hash="none" records only the path.
%
% INPUTS (name-value):
%   Script    the script or function that produced the result; defaults to the caller's file
%   Data      string array of data file paths (rheome.load.root() entries, caches, raw recordings)
%   Extra     further source files to pin (helpers the script depends on)
%   Hash      "fingerprint" (default) | "sha256" | "none"
%   Repo      repository root; defaults to the working directory
%
% OUTPUT (struct P):
%   .git      .branch .sha .shaShort .clean .dirtyFiles .subject .commitDate .remote
%   .code     table: file, role ("script"|"extra"), bytes, modified (formatted string), digest
%   .data     table: file, bytes, modified, digest, exists
%   .env      .matlab .release .os .host .user .toolboxes (table: name, version)
%   .when     datetime the report was generated
%   .table    the whole thing flattened to one [key value] table, for a document
%
% See also: rheome.report.article, rheome.report.template
%
% Author: Diellor Basha, 2026

    arguments
        opts.Script (1,1) string = ""
        opts.Data         string = string.empty
        opts.Extra        string = string.empty
        opts.Hash   (1,1) string {mustBeMember(opts.Hash, ["fingerprint","sha256","none"])} = "fingerprint"
        opts.Repo   (1,1) string = string(pwd)
    end

    if opts.Script == ""
        st = dbstack('-completenames');
        if numel(st) > 1, opts.Script = string(st(2).file); else, opts.Script = "<console>"; end
    end

    P = struct();
    P.when = datetime('now', 'TimeZone', 'local');
    P.git  = i_git(opts.Repo);

    code = [struct('file', opts.Script, 'role', "script")];
    for i = 1:numel(opts.Extra)
        code(end+1) = struct('file', opts.Extra(i), 'role', "extra");  %#ok<AGROW>
    end
    P.code = i_files(string({code.file})', string({code.role})', opts.Hash);
    P.data = i_files(opts.Data(:), repmat("data", numel(opts.Data), 1), opts.Hash);

    v = ver;
    P.env = struct( ...
        'matlab',  string(version), ...
        'release', string(version('-release')), ...
        'os',      string(computer('arch')), ...
        'host',    i_host(), ...
        'user',    string(getenv('USER')), ...
        'toolboxes', table(string({v.Name})', string({v.Version})', ...
                           'VariableNames', {'name','version'}));
    P.table = i_flatten(P);
end

% ---------------------------------------------------------------------------------------
function g = i_git(root)
    g = struct('branch',"", 'sha',"", 'shaShort',"", 'clean',false, 'dirtyFiles',string.empty, ...
               'subject',"", 'commitDate',"", 'remote',"", 'available',false);
    % ⚠ GIT PAGES EVEN THROUGH system(). Without --no-pager the subject and date come back
    % wrapped in terminal escape sequences ("<ESC>[?1h=feat(show): ...<ESC>[m<ESC>[K<ESC>[?1l>")
    % and land in the document looking like corruption. i_run strips whatever survives.
    q = @(c) i_run(sprintf('git --no-pager -C "%s" %s', root, c));
    [ok, ~] = q('rev-parse --is-inside-work-tree');
    if ~ok, return; end
    g.available = true;
    [~, g.branch]     = q('rev-parse --abbrev-ref HEAD');
    [~, g.sha]        = q('rev-parse HEAD');
    [~, g.shaShort]   = q('rev-parse --short HEAD');
    [~, g.subject]    = q('log -1 --pretty=%s');
    [~, g.commitDate] = q('log -1 --pretty=%cI');
    [~, g.remote]     = q('config --get remote.origin.url');
    [~, st]           = q('status --porcelain');
    st = strtrim(st);
    if st == ""
        g.clean = true;
    else
        g.clean = false;
        g.dirtyFiles = strtrim(split(st, newline));
    end
end

function [ok, out] = i_run(cmd)
    [s, o] = system(cmd);
    ok = (s == 0);
    o = regexprep(o, '\x1b\[[0-9;?]*[a-zA-Z=]', '');     % CSI sequences
    o = regexprep(o, '[\x00-\x08\x0b\x0c\x0e-\x1f]', '');
    out = string(strtrim(o));
    if ~ok, out = ""; end
end

function T = i_files(files, roles, hashMode)
    n = numel(files);
    % ⚠ a datetime column renders as "[1x1 datetime]" in a MATLABTable, so the timestamp is
    % formatted here rather than in the document
    T = table(strings(n,1), strings(n,1), zeros(n,1), strings(n,1), strings(n,1), false(n,1), ...
              'VariableNames', {'file','role','bytes','modified','digest','exists'});
    for i = 1:n
        f = files(i);
        T.file(i) = f;  T.role(i) = roles(i);
        d = dir(f);
        if isempty(d) || isfolder(f), continue; end
        T.exists(i)   = true;
        T.bytes(i)    = d.bytes;
        T.modified(i) = string(datetime(d.datenum, 'ConvertFrom', 'datenum', ...
                                        'Format', 'yyyy-MM-dd HH:mm'));
        switch hashMode
            case "sha256",      T.digest(i) = i_sha256(f);
            case "fingerprint", T.digest(i) = sprintf('%d:%s', d.bytes, datestr(d.datenum, 'yyyymmddHHMMSS'));
            otherwise,          T.digest(i) = "";
        end
    end
end

function h = i_sha256(f)
    fid = fopen(f, 'r');  c = onCleanup(@() fclose(fid));
    md = java.security.MessageDigest.getInstance('SHA-256');
    while true
        b = fread(fid, 1e7, '*uint8');       % 10 MB at a time: a 300 MB cache is not read whole
        if isempty(b), break; end
        md.update(b);
    end
    d = typecast(md.digest(), 'uint8');
    h = string(lower(reshape(dec2hex(d, 2)', 1, [])));
end

function h = i_host()
    [ok, o] = i_run('hostname');
    if ok, h = o; else, h = "unknown"; end
end

function T = i_flatten(P)
    k = ["generated"; "git branch"; "git commit"; "git subject"; "committed"; "working tree"; ...
         "MATLAB"; "platform"; "host"];
    dirtyTxt = "clean";
    if P.git.available && ~P.git.clean
        dirtyTxt = sprintf("DIRTY -- %d uncommitted file(s)", numel(P.git.dirtyFiles));
    elseif ~P.git.available
        dirtyTxt = "not a git repository";
    end
    v = [string(P.when, 'yyyy-MM-dd HH:mm:ss zzz'); P.git.branch; P.git.sha; P.git.subject; ...
         P.git.commitDate; dirtyTxt; P.env.release; P.env.os; P.env.host];
    T = table(k, v, 'VariableNames', {'item','value'});
end
% Author: Diellor Basha, 2026
