function rheome_mcp_call(jobDir)
% RHEOME_MCP_CALL  Run one Rheome MCP call: read jobDir/call.json, call its function, write jobDir/result.json.
%
%   matlab -batch "addpath('<mcp>'); rheome_mcp_call('<jobDir>')"      (started by rheome_mcp.py, never by hand)
%
% call.json: root (the Rheome install to put on the path), function (one of the server's whitelist, without the
% "rheome." prefix), args (positional), options (name-value), returns (true: jsonencode the first output into
% result.json). result.json: ok, message, matlab release, resolved (which() of the function -- proof of which
% install ran), seconds, value.

    t0 = tic;
    R = struct('ok', false, 'message', '', 'matlab', version, 'resolved', '', 'seconds', NaN, 'value', []);
    try
        c = jsondecode(fileread(fullfile(jobDir, 'call.json')));
        addpath(c.root);
        fn = "rheome." + string(c.function);
        R.resolved = which(char(fn));
        pos = i_cell(c.args);
        nv = {};
        if isstruct(c.options)
            for f = fieldnames(c.options)'
                nv = [nv, f, {c.options.(f{1})}]; %#ok<AGROW>
            end
        end
        if c.returns
            v = feval(fn, pos{:}, nv{:});
            R.value = v;
        else
            feval(fn, pos{:}, nv{:});
        end
        R.ok = true;
    catch e
        R.message = getReport(e, 'extended', 'hyperlinks', 'off');
    end
    R.seconds = toc(t0);
    fid = fopen(fullfile(jobDir, 'result.json.tmp'), 'w');
    fwrite(fid, jsonencode(R), 'char');
    fclose(fid);
    movefile(fullfile(jobDir, 'result.json.tmp'), fullfile(jobDir, 'result.json'));   % atomic: the server never reads half a file
end

function c = i_cell(a)
% jsondecode gives [] for an empty list, a cell for mixed lists, and an array for a uniform one
    if isempty(a), c = {}; elseif iscell(a), c = a(:)'; else, c = num2cell(a(:)'); end
end
