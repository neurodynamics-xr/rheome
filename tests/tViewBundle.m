classdef tViewBundle < matlab.unittest.TestCase
% TVIEWBUNDLE  rheome.io.write.bundle writes rheome.bundle@1 as the viewer reads it; rheome.view maps its options to the
% viewer's settings and needs no viewer: 'none' writes beneath RHEOME_VIEWER (one bundle at a time), 'show' draws with
% rheome.show. On R2026a+ with a display, the uihtml path is run against a stand-in viewer.html that answers with the
% settings it received (round trip, a refused setting as a warning, figure reuse, set, screenshot, the 'no-webgpu'
% fallback). The real viewer (nxr-cortical-flow editor/dist-viewer) is checked by hand: it needs WebGPU.
%
% Author: Diellor Basha, 2026
    properties
        S = struct('Vertices', [0 0 0; 1 0 0; 0 1 0; 0 0 1], 'Faces', [1 2 3; 1 3 4; 1 4 2; 2 4 3]);
    end
    methods (Test)
        function writesTheFormat(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));
            I = rand(4, 3);  v = rand(12, 3);
            m = rheome.io.write.bundle(d, tc.S, struct('activity', I, 'velocity', v), 'Fs', 100);
            j = jsondecode(fileread(m));
            tc.verifyEqual(j.format, 'rheome.bundle@1');
            tc.verifyEqual([j.time.samples j.time.fs], [3 100]);
            tc.verifyEqual({j.fields.type}, {'real', 'vector'});
            tc.verifyEqual(reshape(rd(d, 'positions.bin', 'single'), 3, []), single(tc.S.Vertices'));
            tc.verifyEqual(reshape(rd(d, 'faces.bin', 'int32'), 3, []), int32(tc.S.Faces' - 1));
            tc.verifyEqual(reshape(rd(d, j.fields(2).values, 'single'), 12, 3), single(v));
        end
        function refusesABadResult(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));  mkdir(d);
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, rand(5, 2)), 'rheome:bundle:rows');
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, struct('a', rand(4, 2), 'b', rand(4, 3))), 'rheome:bundle:time');
        end
        function noneWritesBeneathTheViewer(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            a = rheome.view(tc.S, rand(4, 2), 'Viewer', 'none');
            pause(0.002);
            b = rheome.view(tc.S, rand(12, 2), 'Viewer', 'none');
            tc.verifyTrue(startsWith(b.manifest, fullfile(d, 'bundles')));
            tc.verifyTrue(isfile(b.manifest));
            tc.verifyFalse(isfile(a.manifest));   % one bundle at a time
            tc.verifyEqual(b.viewer, 'none');
        end
        function showDrawsWithoutTheViewer(tc)
            out = rheome.view(tc.S, struct('speed', rand(4, 2)), 'Viewer', 'show');
            c = onCleanup(@() close(out.figure));
            tc.verifyEqual(out.viewer, 'show');
            tc.verifyClass(out.figure, 'matlab.ui.Figure');
            tc.verifyEmpty(out.manifest);
        end
        function writesSurfaces(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));
            Si = tc.S;  Si.Vertices = 2 * tc.S.Vertices;
            m = rheome.io.write.bundle(d, struct('pial', tc.S, 'inflated', Si), rand(4, 2));
            j = jsondecode(fileread(m));
            tc.verifyEqual(j.surface.name, 'pial');
            tc.verifyEqual(j.surfaces(1).name, 'inflated');
            tc.verifyEqual(reshape(rd(d, j.surfaces(1).positions, 'single'), 3, []), single(Si.Vertices'));
            tc.verifyFalse(isfield(j, 'settings'));   % no settings given, none written
        end
        function mapsEveryOption(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            r = struct('activity', rand(4, 5), 'velocity', rand(12, 5));
            out = rheome.view(struct('pial', tc.S, 'inflated', tc.S), r, Viewer="none", Fs=10, ...
                Colormap="viridis", CLim=[-1 2], Threshold=[0.1 NaN], Surface="Inflated", ...
                Hemisphere="left", View="Top", ArrowDensity=0.25, ArrowScale=0.002, Time=0.26, Play=false, Speed=0.5);
            want = struct('layers', {{struct('field', 'activity', 'form', 'f', 'as', 'colour', 'colormap', 'viridis', ...
                'clim', [-1 2], 'threshold', struct('min', 0.1, 'max', NaN)), ...
                struct('field', 'velocity', 'form', 'x', 'as', 'arrows', 'density', 0.25, 'scale', 0.002)}}, ...
                'surface', 'inflated', 'hemisphere', 'left', 'view', 'top', 'frame', 3, 'playing', false, 'fps', 5);
            tc.verifyEqual(out.settings, want);
            txt = fileread(out.manifest);
            tc.verifySubstring(txt, '"threshold":{"min":0.1,"max":null}');   % NaN: the open side
            j = jsondecode(txt);
            tc.verifyEqual(j.settings.layers{1}.clim, [-1; 2]);
            tc.verifyFalse(isfield(j.settings, 'scalar') || isfield(j.settings, 'vector'));   % the old keys are not sent
            o4 = rheome.view(tc.S, rand(4, 2), Viewer="none");
            tc.verifyEqual(o4.settings, struct());   % none given: the viewer's defaults
        end
        function refusesBadOptions(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            v = @(varargin) rheome.view(tc.S, struct('a', rand(4, 3)), 'Viewer', 'none', varargin{:});
            tc.verifyError(@() v('Layers', {{'b'}}), 'rheome:view:option');
            tc.verifyError(@() v('Layers', {{'a', 'abs'}}), 'rheome:view:option');               % a real field has form f only
            tc.verifyError(@() v('Layers', {{'a', 'f', 'arrows'}}), 'rheome:view:option');       % a scalar form is not arrows
            tc.verifyError(@() v('Layers', struct('field', 'a', 'size', 2)), 'rheome:view:option');   % not a layer key
            tc.verifyError(@() v('Layers', struct('field', 'a', 'opacity', 2)), 'rheome:view:option');
            tc.verifyError(@() v('CLim', [1 0]), 'rheome:view:option');
            tc.verifyError(@() v('Colormap', gray(4)), 'rheome:view:option');   % the viewer takes names
            tc.verifyError(@() v('Threshold', [1 0]), 'rheome:view:option');
            tc.verifyError(@() v('ArrowScale', 1), 'rheome:view:option');       % no arrows layer to act on
            tc.verifyError(@() v('Surface', 'inflated'), 'rheome:view:option');
            tc.verifyError(@() v('View', 'lateral'), 'rheome:view:option');
            tc.verifyError(@() v('ArrowDensity', 1.5), 'rheome:view:option');
            tc.verifyError(@() v('Time', 3), 'rheome:view:option');                % T = 3 at 1 Hz: [0, 2] s
            tc.verifyError(@() v('Speed', 0), 'rheome:view:option');
            tc.verifyError(@() rheome.view(rand(4, 3)), 'rheome:view:input');
        end
        function takesTheResultAlone(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            out = rheome.view(struct('surface', tc.S, 'power', rand(4, 2)), Viewer="none", Layers={{"power"}});
            j = jsondecode(fileread(out.manifest));
            tc.verifyEqual({j.fields.name}, {'power'});
            tc.verifyEqual(out.settings.layers, {struct('field', 'power', 'form', 'f', 'as', 'colour')});
        end
        function showHonoursFieldAndTime(tc)
            r = struct('a', zeros(4, 3), 'b', [zeros(4, 2) (1:4)']);
            out = rheome.view(tc.S, r, Viewer="show", Layers=struct('field', 'b', 'clim', [0 4]), Time=2, View="top");
            c = onCleanup(@() close(out.figure));
            ax = findobj(out.figure, 'Type', 'axes');
            p = findobj(ax, 'Type', 'patch');
            tc.verifyEqual(p(1).FaceVertexCData, (1:4)');
            tc.verifyEqual(ax(1).CLim, [0 4]);
            tc.verifyEqual(ax(1).View, [0 90]);
        end
        function takesLayers(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            r = struct('activity', rand(4, 2), 'J', rand(12, 2));
            o = rheome.view(tc.S, r, Viewer="none", Layers={{"J", "curl", "colour"}, {"J", "curl", "halfMax"}, {"j"}});
            tc.verifyEqual(o.settings.layers, {struct('field', 'J', 'form', 'curl', 'as', 'colour'), ...
                struct('field', 'J', 'form', 'curl', 'as', 'halfMax'), struct('field', 'J', 'form', 'x', 'as', 'arrows')});
            o = rheome.view(tc.S, r, Viewer="none", Layers={"J", "abs"});              % one layer; as: colour
            tc.verifyEqual(o.settings.layers, {struct('field', 'J', 'form', 'abs', 'as', 'colour')});
            L = struct('field', {'activity', 'J'}, 'as', {'colour', 'arrows'}, 'color', {[], 'paper'});
            o = rheome.view(tc.S, r, Viewer="none", Layers=L);
            tc.verifyEqual(o.settings.layers{2}.color, 'paper');
            j = jsondecode(fileread(o.manifest));
            tc.verifyEqual(cellfun(@(x) x.field, j.settings.layers, 'UniformOutput', false)', {'activity', 'J'});           % a JSON array, set whole
            tc.verifyError(@() rheome.view(tc.S, r, Viewer="none", Layers={{"J", "x"}}, CLim=[0 1]), 'rheome:view:option');   % no colour layer
        end
        function keepsFieldAndArrowsAsAliases(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            r = struct('activity', rand(4, 2), 'speed', rand(4, 2), 'velocity', rand(12, 2));
            o = tc.verifyWarning(@() rheome.view(tc.S, r, Viewer="none", Field="speed"), 'rheome:view:deprecated');
            tc.verifyEqual(o.settings.layers, {struct('field', 'speed', 'form', 'f', 'as', 'colour'), ...
                struct('field', 'velocity', 'form', 'x', 'as', 'arrows')});
            o = tc.verifyWarning(@() rheome.view(tc.S, r, Viewer="none", Arrows=false), 'rheome:view:deprecated');
            tc.verifyEqual(o.settings.layers, {struct('field', 'activity', 'form', 'f', 'as', 'colour')});
            o = tc.verifyWarning(@() rheome.view(tc.S, r, Viewer="none", Field="velocity"), 'rheome:view:deprecated');
            tc.verifyEqual(o.settings.layers{2}.field, 'velocity');
            tc.verifyError(@() rheome.view(tc.S, struct('a', rand(4, 2)), Viewer="none", Arrows=true), 'rheome:view:option');
        end
        function writesPotentials(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));
            [V, F] = rheome.geom.icosphere(2);  S = struct('Vertices', 0.07*V, 'Faces', F);
            nV = size(V, 1);  J = randn(3*nV, 2);
            m = rheome.io.write.bundle(d, S, struct('a', rand(nV, 2), 'J', J), 'Potentials', "J");
            j = jsondecode(fileread(m));
            tc.verifyEqual(cellfun(@(x) x.name, j.fields, 'UniformOutput', false)', {'a', 'J', 'J_phi', 'J_psi'});
            g = j.fields(3:4);  if iscell(g), g = [g{:}]; end
            tc.verifyEqual({g.name; g.derivedFrom; g.form; g.type}, {'J_phi', 'J_psi'; 'J', 'J'; 'phi', 'psi'; 'real', 'real'});
            fg = rheome.operators.face_gradient(S.Vertices, S.Faces);  S.VertNormals = fg.W * fg.FaceNormal;
            H = rheome.differential.helmholtz(J, S);
            tc.verifyEqual(reshape(rd(d, g(2).values, 'single'), nV, 2), single(H.Psi));
            tc.verifyError(@() rheome.io.write.bundle(d, S, struct('J', J), 'Potentials', "a"), 'rheome:bundle:potentials');
            old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', [d '_v']);  c2 = onCleanup(@() cleanup([d '_v'], old));
            out = rheome.view(S, struct('J', J), Viewer="none", Potentials=true, Layers={{"J_phi"}, {"J"}});
            tc.verifyEqual(out.settings.layers{1}, struct('field', 'J_phi', 'form', 'f', 'as', 'colour'));
            tc.verifyTrue(contains(fileread(out.manifest), '"derivedFrom":"J"'));
        end
        function editorRoundTrip(tc)
            tc.assumeFalse(isMATLABReleaseOlderThan('R2026a'), 'uihtml viewer path: R2026a+');
            tc.assumeTrue(usejava('desktop') || feature('ShowFigureWindows'), 'needs a display');
            d = tempname;  mkdir(d);  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            writePage(d, 'echo');
            a = tc.verifyWarning(@() rheome.view(tc.S, struct('a', rand(4, 3)), Viewer="editor", CLim=[0 1], Time=1, ...
                Timeout=30), 'rheome:view:settings');   % the page refuses one setting first: a warning, then ready
            tc.addTeardown(@() delete(a.figure));
            tc.verifyEqual(a.viewer, 'editor');
            tc.verifyEqual(a.ready.settings.layers.clim(:)', [0 1]);
            tc.verifyEqual(a.ready.settings.frame, 1);
            tc.verifyEqual(a.ready.bytes(:)', [48 48 48]);   % positions, faces, the field: every file served
            b = rheome.view(tc.S, struct('a', rand(4, 3)), Viewer="editor", Layers={{"a", "f", "isolines"}}, Fs=100, Timeout=30);
            tc.verifySameHandle(b.figure, a.figure);    % the figure is reused
            tc.verifyEqual(b.ready.settings.layers.as, 'isolines');
            h = findobj(b.figure, 'Type', 'uihtml');
            b.set(View="right", Play=true, Speed=0.5);
            t = tic;  while ~isfield(h.UserData, 'settings') && toc(t) < 30, pause(0.1); end
            tc.verifyEqual(h.UserData.settings, struct('view', 'right', 'playing', true, 'fps', 50));   % Speed x the bundle's Fs
            tc.verifyError(@() b.set(View="sideways"), 'rheome:view:option');
            h.UserData = rmfield(h.UserData, 'settings');
            b.set(CLim=[0 2]);   % a layer option acts on the layers last sent
            t = tic;  while ~isfield(h.UserData, 'settings') && toc(t) < 30, pause(0.1); end
            tc.verifyEqual(h.UserData.settings.layers, struct('field', 'a', 'form', 'f', 'as', 'isolines', 'clim', [0 2]));
            f = [tempname '.png'];  tc.addTeardown(@() delete(f));
            b.screenshot(f, 2);
            tc.verifyEqual(size(imread(f)), [1 1 3]);
        end
        function editorWithoutWebGPUShows(tc)
            tc.assumeFalse(isMATLABReleaseOlderThan('R2026a'), 'uihtml viewer path: R2026a+');
            tc.assumeTrue(usejava('desktop') || feature('ShowFigureWindows'), 'needs a display');
            d = tempname;  mkdir(d);  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            writePage(d, 'nogpu');
            out = rheome.view(tc.S, rand(4, 2), Timeout=30);
            tc.addTeardown(@() close(out.figure));
            tc.verifyEqual(out.viewer, 'show');
            tc.verifyClass(out.figure, 'matlab.ui.Figure');
        end
    end
end

function writePage(d, kind)
% a stand-in for the viewer's viewer.html with its protocol: Data {bundle, settings}, DataChanged, set, screenshot;
% events ready, settings, screenshot, error (one refused setting before each ready, as the viewer does)
    if strcmp(kind, 'nogpu')
        js = 'function setup(h){ h.sendEventToMATLAB(''error'', {code: ''no-webgpu'', message: ''no adapter''}) }';
    else
        js = ['function setup(h){ const go = async () => { try {' ...
            ' const u = new URL(h.Data.bundle, location.href), m = await (await fetch(u)).json();' ...
            ' const files = [m.surface.positions, m.surface.faces, ...m.fields.map((f) => f.values)];' ...
            ' const bytes = await Promise.all(files.map(async (f) => { const r = await fetch(new URL(f, u));' ...
            ' if (!r.ok) throw new Error(f + '': '' + r.status); return (await r.arrayBuffer()).byteLength }));' ...
            ' h.sendEventToMATLAB(''error'', {code: ''settings'', message: ''bogus: not a setting''});' ...
            ' h.sendEventToMATLAB(''ready'', {settings: h.Data.settings, bytes}) }' ...
            ' catch (e) { h.sendEventToMATLAB(''error'', {code: ''bundle'', message: String(e)}) } };' ...
            ' h.addEventListener(''set'', (e) => h.sendEventToMATLAB(''settings'', e.Data));' ...
            ' h.addEventListener(''screenshot'', () => h.sendEventToMATLAB(''screenshot'', {width: 1, height: 1, png:' ...
            ' ''data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAIAAACQd1PeAAAADElEQVR4nGP4z8AAAAMBAQDJ/pLvAAAAAElFTkSuQmCC''}));' ...
            ' h.addEventListener(''DataChanged'', go); go() }'];
    end
    fid = fopen(fullfile(d, 'viewer.html'), 'w');
    fprintf(fid, '<!doctype html><html><body><script>%s</script></body></html>', js);
    fclose(fid);
end

function x = rd(d, f, cls)
    fid = fopen(fullfile(d, f), 'r', 'ieee-le');  x = fread(fid, inf, ['*' cls]);  fclose(fid);
end

function cleanup(d, old)
    setenv('RHEOME_VIEWER', old);
    if isfolder(d), rmdir(d, 's'); end
end
