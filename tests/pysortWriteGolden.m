function pysortWriteGolden(txt)
% Converts pysortCaptureGolden.py's text (one case per line: n, n keys, n 1-based sorted indices)
% into tests/pysortGolden.mat for tSurfaceLadder.
lines = splitlines(strtrim(fileread(txt)));
keys = cell(numel(lines), 1);  perm = keys;
for i = 1:numel(lines)
    v = str2double(split(lines{i}));  n = v(1);
    keys{i} = v(2:n+1);  perm{i} = v(n+2:end);
end
save(fullfile(fileparts(mfilename('fullpath')), 'pysortGolden.mat'), 'keys', 'perm');
end
