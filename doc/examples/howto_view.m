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

%% Choose what is shown: layers
% The painted colours and the arrows are both fields, of two types: a real
% field (a number per vertex) and a vector field. What is drawn is a list
% of |Layers|, in draw order. Each layer is one field, one *form* of it and
% one *representation*: a real field's value |"f"|; a vector field's
% strength |"abs"|, normal part |"nx"|, divergence |"div"| or rotation
% |"curl"| (scalar forms, drawn as |"colour"|, |"isolines"| or
% |"halfMax"|), or the vector itself |"x"| (drawn as |"arrows"|). A vector
% field is never painted itself. Here the bump's activity is painted, its
% half-maximum outlined, and its flow drawn as arrows:

out = rheome.view(S, result, Fs=fs, Layers={{"activity", "f", "colour"}, ...
    {"activity", "f", "halfMax"}, {"velocity", "x", "arrows"}});
out.settings.layers{:}

%%
% The rotation of the flow, painted, with the flow over it:

out = rheome.view(S, result, Fs=fs, Layers={{"velocity", "curl"}, {"velocity"}});

%%
% A layer's look is set in a struct with the viewer's layer keys
% (|colormap|, |clim|, |threshold|, |opacity|, |density|, |color|, ...).
% |Colormap|, |CLim| and |Threshold| set the colour layer's look, and
% |ArrowDensity| and |ArrowScale| every arrows layer's. The other options
% are the surface, the camera and the time. An option you do not give is
% left to the viewer. A wrong one (a field that is not in the result, a
% form the field does not have, |CLim| with lo >= hi, a time after the last
% sample) is an error before anything is written. |Threshold| is a band:
% only values in [min, max] are painted (NaN leaves a side open).

out = rheome.view(S, result, Fs=fs, ...
    Layers=struct('field', {'activity', 'velocity'}, 'as', {'colour', 'arrows'}, 'color', {[], 'paper'}), ...
    Colormap="viridis", CLim=[0 1], Threshold=[0.05 NaN], ...
    View="top", ArrowDensity=0.2, Time=0.1, Play=true, Speed=0.25);

%% Sources and vortices: the Hodge potentials
% The scalar potential phi (its maxima are sources, its minima sinks) and
% the stream function psi (its extrema are vortices) of a vector field
% need a Poisson solve, which the viewer does not do. |Potentials| has
% Rheome solve them (|rheome.differential.helmholtz|) and write them as
% real fields |<v>_phi| and |<v>_psi|, tagged with the field they come
% from:

out = rheome.view(S, result, Fs=fs, Potentials="velocity", ...
    Layers={{"velocity_psi"}, {"velocity", "x", "arrows"}});

%% Change the open view, and save it
% Inside MATLAB, |out.set| changes the open view with the same options,
% without reopening it, and |out.screenshot| writes it as a PNG:
%
%   out.set(View="left", Play=false, Time=0.05)
%   out.screenshot("bump.png", 2)   % 2: twice the window's pixels

%% Several surfaces, one hemisphere
% Pass a struct of surfaces of one mesh, by name (|pial|, |inflated|,
% |sphere|, ...), and pick one with |Surface|. |Hemisphere| shows one of
% the mesh's two largest connected components (the one at larger y is the
% left, as in Brainstorm's SCS); this sphere is one piece, so it is not
% used here:
%
%   surfs = struct('pial', P, 'inflated', Pi);
%   rheome.view(surfs, result, Surface="inflated", Hemisphere="left", View="left")

surfs = struct('pial', S, 'sphere', setfield(S, 'Vertices', 0.1*V)); %#ok<SFLD>
out = rheome.view(surfs, result, Fs=fs, Surface="sphere", View="front");

%% Reuse the window
% Inside MATLAB the next call opens its result in the same viewer window,
% in place of the last one; close the window to get a new one. Only the
% last bundle is kept on disk (beneath the viewer's folder).
%
%% The result alone
% A struct that carries its surface in a |surface| field is passed alone:

r = result;  r.surface = S;
out = rheome.view(r, Fs=fs, Layers={{"velocity", "abs"}, {"velocity"}});

%% Field and Arrows (deprecated)
% |Field| and |Arrows|, the options before layers, still work for this
% release, with a warning, and are removed in the next. |Field="x"| is a
% layer |{"x"}|; |Arrows=false| drops the arrows layers.

%% Install the viewer
% Inside MATLAB (R2026a+) the viewer is a folder with |viewer.html| (the
% Cortical Flow viewer's build, |editor/dist-viewer|).
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
