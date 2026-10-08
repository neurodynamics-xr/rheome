function check_install(mltbx)
% CHECK_INSTALL  Install the .mltbx into a clean path, run the quick start, uninstall.
%
%   check_install('build/out/rheome-1.0.0.mltbx')
%
% Run it in a FRESH MATLAB started outside the repository, so the clone is not on the path:
%   cd /tmp && matlab -batch "addpath('<repo>/build'); check_install('<repo>/build/out/<name>.mltbx')"
% It restores the default path (build/ goes too), installs the add-on, and asserts:
% every name resolves inside the add-on folder; the quick start (README) passes its self-checks;
% rheome.load.root/outroot default under userpath, not inside the add-on. Then it uninstalls.
%
% See also: build_toolbox
%
% Author: Diellor Basha, 2026

    mltbx = char(java.io.File(mltbx).getCanonicalPath());   % absolute, before the path changes
    restoredefaultpath;
    assert(isempty(which('rheome.graphfilterbank')), 'check_install:dirty', ...
        'rheome.graphfilterbank is already on the path (%s): start MATLAB outside the clone.', which('rheome.graphfilterbank'));
    setenv('RHEOME_DATA', '');  setenv('RHEOME_OUT', '');

    tb = matlab.addons.install(mltbx);
    cleanup = onCleanup(@() matlab.addons.uninstall(tb.Identifier));
    fprintf('installed %s %s\n', tb.Name, tb.Version);

    addon = fileparts(fileparts(fileparts(which('rheome.graphfilterbank'))));   % <addon>/+rheome/@graphfilterbank/x.m
    fprintf('add-on folder: %s\n', addon);
    for n = ["rheome.load.root", "rheome.flow.context", "rheome.demos.filterbank_graph", "rheome.jointfilterbank"]
        assert(startsWith(which(n), addon), 'check_install:resolve', '%s resolves to %s', n, which(n));
    end

    % --- the README quick start, verbatim ---
    out = rheome.demos.filterbank_graph;
    assert(out.ok, 'check_install:quickstart', 'rheome.demos.filterbank_graph self-checks failed');
    [V, F] = rheome.geom.icosphere(4);
    [L, M] = rheome.operators.laplace_beltrami(V, F, 'galerkin');
    gfb = rheome.graphfilterbank.fromOperator(L, 'Mass', M);
    x = double(V(:,3) > 0.9);
    E = vertexSpectrum(gfb, x);
    assert(numel(E) == gfb.NumMembers && all(isfinite(E)) && any(E > 0), 'check_install:quickstart', ...
        'vertexSpectrum of the polar cap is not a finite, non-zero [M x 1] spectrum');
    disp(E.')
    close all

    % --- outputs never default into the add-on folder ---
    u = fullfile(userpath, 'rheome');
    assert(strcmp(rheome.load.root(), fullfile(u, 'data')), 'check_install:root', 'rheome.load.root = %s', rheome.load.root());
    assert(strcmp(rheome.load.outroot(), fullfile(u, 'results')), 'check_install:outroot', 'rheome.load.outroot = %s', rheome.load.outroot());
    fprintf('rheome.load.root    = %s\nrheome.load.outroot = %s\n', rheome.load.root(), rheome.load.outroot());
    fprintf('CHECK_INSTALL PASS\n');
end

% Author: Diellor Basha, 2026
