classdef tViewBundle < matlab.unittest.TestCase
% TVIEWBUNDLE  rheome.io.write.bundle writes rheome.bundle@1 as the viewer reads it; rheome.view maps its options to the
% viewer's settings and needs no viewer: 'none' writes beneath RHEOME_VIEWER (one bundle at a time), 'show' draws with
% rheome.show. On R2026a+ with a display, the uihtml path is run against a stand-in viewer page that answers with the
% settings it received (round trip, figure reuse, the 'no-webgpu' fallback).
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
        function writesSurfacesAndHemispheres(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));
            S = tc.S;  S.Hemi = {[1 2]', [3 4]'};  S.HemiLabel = {'Cortex L', 'Cortex R'};
            Si = S;  Si.Vertices = 2 * S.Vertices;
            m = rheome.io.write.bundle(d, struct('pial', S, 'inflated', Si), rand(4, 2));
            j = jsondecode(fileread(m));
            tc.verifyEqual(j.surface.name, 'pial');
            tc.verifyEqual(j.surface.alternates(1).name, 'inflated');
            tc.verifyEqual(reshape(rd(d, j.surface.alternates(1).positions, 'single'), 3, []), single(Si.Vertices'));
            tc.verifyEqual(rd(d, j.surface.hemisphere, 'uint8'), uint8([0 0 1 1]'));
            tc.verifyFalse(isfield(j, 'view'));   % no settings given, none written
        end
        function mapsEveryOption(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            S = tc.S;  S.Hemi = {[1 2]', [3 4]'};
            r = struct('activity', rand(4, 5), 'velocity', rand(12, 5), 'critical', rand(4, 5));
            out = rheome.view(struct('pial', S, 'inflated', S), r, Viewer="none", Fs=10, Field="velocity", ...
                Colormap="viridis", CLim=[-1 2], Threshold=0.1, Surface="Inflated", Hemisphere="left", View="medial", ...
                Arrows=true, ArrowDensity=0.25, ArrowScale=2, Overlays=["critical" "activity"], Time=0.3, Play=false, Speed=0.5);
            want = struct('field', 'velocity', 'colormap', 'viridis', 'clim', [-1 2], 'threshold', 0.1, ...
                'surface', 'inflated', 'hemisphere', 'left', 'camera', 'medial', ...
                'arrows', struct('visible', true, 'density', 0.25, 'scale', 2), 'overlays', {{'critical', 'activity'}}, ...
                'time', 0.3, 'play', false, 'speed', 0.5);
            tc.verifyEqual(out.settings, want);
            j = jsondecode(fileread(out.manifest));
            tc.verifyEqual(j.view.clim, [-1; 2]);
            tc.verifyEqual(j.view.overlays, {'critical'; 'activity'});
            tc.verifyEqual(j.view.arrows, want.arrows);
            tc.verifyEqual(rmfield(j.view, {'clim', 'overlays'}), rmfield(want, {'clim', 'overlays'}));
            o2 = rheome.view(tc.S, rand(4, 2), Viewer="none", View=[30 -10], Colormap=gray(4));
            tc.verifyEqual(o2.settings.camera, [30 -10]);
            tc.verifyEqual(jsondecode(fileread(o2.manifest)).view.colormap, gray(4));
            o3 = rheome.view(tc.S, rand(4, 2), Viewer="none");
            tc.verifyEqual(o3.settings, struct());   % none given: the viewer's defaults
        end
        function refusesBadOptions(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            v = @(varargin) rheome.view(tc.S, struct('a', rand(4, 3)), 'Viewer', 'none', varargin{:});
            tc.verifyError(@() v('Field', 'b'), 'rheome:view:option');
            tc.verifyError(@() v('CLim', [1 0]), 'rheome:view:option');
            tc.verifyError(@() v('Colormap', [2 0 0]), 'rheome:view:option');
            tc.verifyError(@() v('Hemisphere', 'left'), 'rheome:view:option');     % no .Hemi
            tc.verifyError(@() v('Surface', 'inflated'), 'rheome:view:option');
            tc.verifyError(@() v('View', 'sideways'), 'rheome:view:option');
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
            tc.verifyEqual(out.settings.field, 'power');
        end
        function showHonoursFieldAndTime(tc)
            r = struct('a', zeros(4, 3), 'b', [zeros(4, 2) (1:4)']);
            out = rheome.view(tc.S, r, Viewer="show", Field="b", Time=2, CLim=[0 4], View="dorsal");
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
            a = rheome.view(tc.S, struct('a', rand(4, 3)), Viewer="editor", CLim=[0 1], Time=1, Timeout=30);
            tc.addTeardown(@() delete(a.figure));
            tc.verifyEqual(a.viewer, 'editor');
            h = findobj(a.figure, 'Type', 'uihtml');
            got = h.UserData.data;
            tc.verifyEqual(got.settings.clim(:)', [0 1]);
            tc.verifyEqual(got.settings.time, 1);
            tc.verifyEqual(got.bytes(:)', [48 48 48]);   % positions, faces, the field: every file served
            b = rheome.view(tc.S, struct('a', rand(4, 3)), Viewer="editor", Field="a", Timeout=30);
            tc.verifySameHandle(b.figure, a.figure);    % the figure is reused
            tc.verifyEqual(h.UserData.data.settings.field, 'a');
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
% a stand-in for the viewer's matlab.html with its protocol: Data {bundle, settings}, DataChanged, events ready/error
    if strcmp(kind, 'nogpu')
        js = 'function setup(h){ h.sendEventToMATLAB(''error'', {code: ''no-webgpu'', message: ''no adapter''}) }';
    else
        js = ['function setup(h){ const go = async () => { try {' ...
            ' const u = new URL(h.Data.bundle, location.href), m = await (await fetch(u)).json();' ...
            ' const files = [m.surface.positions, m.surface.faces, ...m.fields.map((f) => f.values)];' ...
            ' const bytes = await Promise.all(files.map(async (f) => { const r = await fetch(new URL(f, u));' ...
            ' if (!r.ok) throw new Error(f + '': '' + r.status); return (await r.arrayBuffer()).byteLength }));' ...
            ' h.sendEventToMATLAB(''ready'', {settings: h.Data.settings, bytes}) }' ...
            ' catch (e) { h.sendEventToMATLAB(''error'', {message: String(e)}) } };' ...
            ' h.addEventListener(''DataChanged'', go); go() }'];
    end
    fid = fopen(fullfile(d, 'matlab.html'), 'w');
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
