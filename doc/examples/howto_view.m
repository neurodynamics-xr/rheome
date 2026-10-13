%% How to view a result in 3D
% *How-to guide.* A travelling pattern is easier to judge moving, on the
% cortex, than frame by frame. |rheome.view| opens a result (a field on the
% cortex over time, and its flow) in Rheome's 3D viewer in one call: you can
% play it, turn the brain, and change the colours. The viewer is *optional*.
% It needs MATLAB R2026a or later (inside MATLAB) or the Cortical Flow
% desktop app (on older releases). Without either, or without WebGPU,
% |rheome.view| draws one frame with |rheome.show|, and nothing else in
% Rheome needs the viewer.

%% A result to view
% This guide uses a synthetic result that runs anywhere: a Gaussian bump
% that moves over a cortex-sized sphere at 100 Hz, and its optical flow.
% Any [nV x T] field, or [3nV x T] vector field, on a surface is viewed
% the same way.

[V, F] = rheome.geom.icosphere(4);
S = struct('Vertices', 0.07*V, 'Faces', F, 'VertNormals', V, 'nV', size(V, 1));
fs = 100;  T = 21;
I = zeros(size(V, 1), T);
for t = 1:T
    c = 0.07 * [cos(0.05*t) sin(0.05*t) 0.3];   % the bump's centre moves along a latitude
    I(:, t) = exp(-sum((S.Vertices - c).^2, 2) / (2*0.012^2));
end
v = rheome.dynamics.opticalflow_scalar(I, S, struct('alpha', 1));   % [3nV x (T-1)]
result = struct('activity', I(:, 1:T-1), 'velocity', v);

%% Open it
% One call writes the result as a bundle (a folder of raw arrays and a
% |manifest.json|, see |rheome.io.write.bundle|) and opens it. |out.viewer|
% says where it opened: |'editor'| (inside MATLAB), |'app'|, or |'show'|
% (the static fallback, as on the machine that built this page).

out = rheome.view(S, result, Fs=fs, Name="Moving bump");
disp(out.viewer)

%% Choose what is shown
% The options are the viewer's settings: the field and the colours, the
% surface and the camera, the arrows, the overlays and the time. An option
% you do not give is left to the viewer. A wrong one (a field that is not
% in the result, |CLim| with lo >= hi, a time after the last sample) is an
% error before anything is written.

out = rheome.view(S, result, Fs=fs, ...
    Field="activity", Colormap="viridis", CLim=[0 1], Threshold=0.05, ...
    View="dorsal", Overlays="velocity", Arrows=true, ArrowDensity=0.2, ArrowScale=1.5, ...
    Time=0.1, Play=true, Speed=0.25);
out.settings

%% Several surfaces, one hemisphere
% Pass a struct of surfaces of one mesh, by name (|pial|, |inflated|,
% |sphere|, ...), and pick one with |Surface|. |Hemisphere| needs the
% hemisphere split that |rheome.io.read.surface| reads from a Brainstorm
% cortex (|.Hemi|); this sphere has none, so it is not used here:
%
%   surfs = struct('pial', P, 'inflated', Pi);
%   rheome.view(surfs, result, Surface="inflated", Hemisphere="left", View="lateral")

surfs = struct('pial', S, 'sphere', setfield(S, 'Vertices', 0.1*V)); %#ok<SFLD>
out = rheome.view(surfs, result, Fs=fs, Surface="sphere", View=[-60 20]);

%% Reuse the window
% Inside MATLAB the next call opens its result in the same viewer window,
% in place of the last one; close the window to get a new one. Only the
% last bundle is kept on disk (beneath the viewer's folder).
%
%% The result alone
% A struct that carries its surface in a |surface| field is passed alone:

r = result;  r.surface = S;
out = rheome.view(r, Fs=fs, Field="velocity", Arrows=true);

%% Install the viewer
% Inside MATLAB (R2026a+) the viewer is a folder with |matlab.html|.
% Rheome uses |RHEOME_VIEWER| if you set it, else |prefdir/rheome-viewer|,
% which it installs from the copy shipped with Rheome (its |viewer| folder)
% when there is one. On R2023b-R2025b, install the Cortical Flow desktop app
% that registers |corticalflow://|. To write the bundle only, for another
% program, pass |Viewer="none"|; to force the static figure, |Viewer="show"|.
%
%% See also
% <reference/rheome.view.html rheome.view>,
% <reference/rheome.io.write.bundle.html rheome.io.write.bundle>,
% <howto_figures.html How to make figures>.
%
% _Written for Rheome @COMMIT@._
