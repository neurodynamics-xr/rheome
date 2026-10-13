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
            out = rheome.view(struct('pial', tc.S, 'inflated', tc.S), r, Viewer="none", Fs=10, Field="activity", ...
                Arrows=false, Colormap="viridis", CLim=[-1 2], Threshold=[0.1 NaN], Surface="Inflated", ...
                Hemisphere="left", View="Top", ArrowDensity=0.25, ArrowScale=0.002, Time=0.26, Play=false, Speed=0.5);
            want = struct('scalar', 'activity', 'vector', '', 'colormap', 'viridis', 'clim', [-1 2], ...
                'threshold', struct('min', 0.1, 'max', NaN), 'surface', 'inflated', 'hemisphere', 'left', 'view', 'top', ...
                'arrowDensity', 0.25, 'arrowScale', 0.002, 'frame', 3, 'playing', false, 'fps', 5);
            tc.verifyEqual(out.settings, want);
            txt = fileread(out.manifest);
            tc.verifySubstring(txt, '"threshold":{"min":0.1,"max":null}');   % NaN: the open side
            j = jsondecode(txt);
            tc.verifyEqual(j.settings.clim, [-1; 2]);
            tc.verifyEqual(j.settings.vector, '');
            o2 = rheome.view(tc.S, r, Viewer="none", Field="velocity");
            tc.verifyEqual(o2.settings, struct('vector', 'velocity'));    % a vector field: the arrows
            o3 = rheome.view(tc.S, r, Viewer="none", Arrows=true);
            tc.verifyEqual(o3.settings, struct('vector', 'velocity'));
            o4 = rheome.view(tc.S, rand(4, 2), Viewer="none");
            tc.verifyEqual(o4.settings, struct());   % none given: the viewer's defaults
        end
        function refusesBadOptions(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            v = @(varargin) rheome.view(tc.S, struct('a', rand(4, 3)), 'Viewer', 'none', varargin{:});
            tc.verifyError(@() v('Field', 'b'), 'rheome:view:option');
            tc.verifyError(@() v('CLim', [1 0]), 'rheome:view:option');
            tc.verifyError(@() v('Colormap', gray(4)), 'rheome:view:option');   % the viewer takes names
            tc.verifyError(@() v('Threshold', [1 0]), 'rheome:view:option');
            tc.verifyError(@() v('Arrows', true), 'rheome:view:option');         % no vector field
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
            out = rheome.view(struct('surface', tc.S, 'power', rand(4, 2)), Viewer="none", Field="power");
            j = jsondecode(fileread(out.manifest));
            tc.verifyEqual({j.fields.name}, {'power'});
            tc.verifyEqual(out.settings.scalar, 'power');
        end
        function showHonoursFieldAndTime(tc)
            r = struct('a', zeros(4, 3), 'b', [zeros(4, 2) (1:4)']);
            out = rheome.view(tc.S, r, Viewer="show", Field="b", Time=2, CLim=[0 4], View="top");
            c = onCleanup(@() close(out.figure));
            ax = findobj(out.figure, 'Type', 'axes');
            p = findobj(ax, 'Type', 'patch');
            tc.verifyEqual(p(1).FaceVertexCData, (1:4)');
            tc.verifyEqual(ax(1).CLim, [0 4]);
            tc.verifyEqual(ax(1).View, [0 90]);
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
            tc.verifyEqual(a.ready.settings.clim(:)', [0 1]);
            tc.verifyEqual(a.ready.settings.frame, 1);
            tc.verifyEqual(a.ready.bytes(:)', [48 48 48]);   % positions, faces, the field: every file served
            b = rheome.view(tc.S, struct('a', rand(4, 3)), Viewer="editor", Field="a", Fs=100, Timeout=30);
            tc.verifySameHandle(b.figure, a.figure);    % the figure is reused
            tc.verifyEqual(b.ready.settings.scalar, 'a');
            h = findobj(b.figure, 'Type', 'uihtml');
            b.set(View="right", Play=true, Speed=0.5);
            t = tic;  while ~isfield(h.UserData, 'settings') && toc(t) < 30, pause(0.1); end
            tc.verifyEqual(h.UserData.settings, struct('view', 'right', 'playing', true, 'fps', 50));   % Speed x the bundle's Fs
            tc.verifyError(@() b.set(View="sideways"), 'rheome:view:option');
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
