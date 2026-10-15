function [file, sha, source] = tokenspath(file, vendored)
% TOKENSPATH  Which figure-tokens file a figure tool reads, and its sha256.
%
%   [file, sha, source] = rheome.viz.tokenspath()
%   [file, sha, source] = rheome.viz.tokenspath(file, vendored)
%
% The order: an explicit FILE; else the environment variable NXR_FIGURE_TOKENS; else the library's
% $HOME/workspace/library/styles/tokens/figure-tokens.json when it exists; else the file vendored with the viewer
% (viewer/figure/figure-tokens.json, or VENDORED). An explicit file or variable that names no file is an error, not
% a silent fall-back. SHA is the sha256 (hex) of the bytes of the file used; SOURCE is "argument" | "env" |
% "library" | "vendored".
%
% Author: Diellor Basha, 2026
    if nargin < 1, file = ''; end
    if nargin < 2 || isempty(vendored)
        vendored = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'viewer', 'figure', 'figure-tokens.json');
    end
    file = char(file);  source = "argument";
    if isempty(file), file = getenv('NXR_FIGURE_TOKENS'); source = "env"; end
    if isempty(file)
        lib = fullfile(getenv('HOME'), 'workspace', 'library', 'styles', 'tokens', 'figure-tokens.json');
        if isfile(lib), file = lib; source = "library"; else, file = char(vendored); source = "vendored"; end
    end
    if ~isfile(file), error('rheome:viz:tokens', 'No figure tokens at %s (Tokens= or NXR_FIGURE_TOKENS).', file); end
    file = char(java.io.File(file).getCanonicalPath());
    fid = fopen(file, 'r');  b = fread(fid, Inf, '*uint8');  fclose(fid);
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(b);
    sha = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, []));
end
