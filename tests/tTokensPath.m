classdef tTokensPath < matlab.unittest.TestCase
% rheome.viz.tokenspath: explicit file > $NXR_FIGURE_TOKENS > library file > vendored copy; sha256 of the file used.
    methods (Test)
        function order(t)
            d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
            a = fullfile(d, 'a.json');  b = fullfile(d, 'b.json');  v = fullfile(d, 'v.json');
            files = {a, '{"x":"a"}'; b, '{"x":"b"}'; v, '{"x":"v"}'};
            for k = 1:3, fid = fopen(files{k,1}, 'w'); fwrite(fid, files{k,2}); fclose(fid); end
            oldE = getenv('NXR_FIGURE_TOKENS');  oldH = getenv('HOME');
            r = onCleanup(@() restore(oldE, oldH));
            setenv('HOME', d);  setenv('NXR_FIGURE_TOKENS', '');
            [f, ~, s] = rheome.viz.tokenspath('', v);
            t.verifyEqual(s, "vendored");  t.verifyEqual(f, char(java.io.File(v).getCanonicalPath()));
            lib = fullfile(d, 'workspace', 'library', 'styles', 'tokens');  mkdir(lib);
            copyfile(a, fullfile(lib, 'figure-tokens.json'));
            [~, ~, s] = rheome.viz.tokenspath('', v);  t.verifyEqual(s, "library");
            setenv('NXR_FIGURE_TOKENS', b);
            [f, ~, s] = rheome.viz.tokenspath('', v);  t.verifyEqual(s, "env");  t.verifyEqual(f, char(java.io.File(b).getCanonicalPath()));
            [f, sha, s] = rheome.viz.tokenspath(v, v);  t.verifyEqual(s, "argument");
            if isunix, [~, w] = system(sprintf('shasum -a 256 "%s"', f)); t.verifyEqual(sha, w(1:64)); end
            t.verifyError(@() rheome.viz.tokenspath(fullfile(d, 'none.json'), v), 'rheome:viz:tokens');
            setenv('NXR_FIGURE_TOKENS', fullfile(d, 'missing.json'));   % a named env file that is missing is an error
            t.verifyError(@() rheome.viz.tokenspath('', v), 'rheome:viz:tokens');
            setenv('NXR_FIGURE_TOKENS', '');  setenv('HOME', '');   % HOME unset: USERPROFILE / java user.home, still resolves
            [~, ~, s] = rheome.viz.tokenspath('', v);  t.verifyTrue(ismember(s, ["library" "vendored"]));
        end
    end
end
function restore(e, h), setenv('NXR_FIGURE_TOKENS', e); setenv('HOME', h); end
