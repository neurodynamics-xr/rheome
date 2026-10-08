%% Getting Started with Rheome
% *Tutorial.* In this tutorial you go from a MEG recording to four maps on the
% cortex: the *divergence* of the cortical current (where it springs and sinks),
% its *rotation* (where it turns), and the two Helmholtz-Hodge potentials,
% the scalar potential Phi and the stream function Psi. You will use a small
% synthetic study that ships with this documentation, so every step runs on
% your machine as written, and at the end you will check that the stream
% function finds the vortex planted in the data.
%
% It takes about a minute. For why each step is what it is, read
% <about_sensor_to_cortex.html From sensors to cortex>; for doing the same
% with your own Brainstorm study, read <howto_ingest_brainstorm.html How to
% import a Brainstorm study>.
%
%% Install the toolbox
% You need MATLAB R2023b or later, with the Signal Processing Toolbox and the
% Statistics and Machine Learning Toolbox (the flow kernels use both).
% Brainstorm itself is not needed: its files are read from disk.
%
% Install the toolbox file by double-clicking |rheome-1.0.0.mltbx|,
% or from the command line:
%
%   matlab.addons.install('rheome-1.0.0.mltbx');
%
% If you work from a clone of the repository instead, put its root on the path:
%
%   addpath('/path/to/rheome');
%
% Check that MATLAB finds the toolbox (2 means a function on the path):

exist('rheome.flow.context')

%% Write the demo study
% The demo study is a Brainstorm study in miniature: a cortex of two
% hemispheres, 150 MEG sensors, an unconstrained leadfield, a 60 s recording
% and a noise covariance. The recording holds a 10 Hz vortex on the lateral
% face of the left hemisphere. The function that writes it lives beside this
% documentation:

root = fileparts(fileparts(fileparts(which('rheome.flow.context'))));   % +rheome/+flow -> toolbox root
addpath(fullfile(root, 'doc', 'examples'))
demo = cfdemo_study();

%% Import it
% Three calls read the study into the toolbox's cache: the cortex and its
% Laplace-Beltrami operators, the channels, leadfield, recording and noise
% covariance, and 100 Laplace-Beltrami eigenmodes per hemisphere.

rheome.import.surface('cfdemo', demo.cortexFile);
rheome.import.study('cfdemo', demo.studyDir, demo.dataName);
rheome.import.bases('cfdemo', 100, 100, struct('overwrite', true));

%% Build the flow context
% |rheome.flow.context| computes the minimum-norm imaging kernel from the leadfield
% and the noise covariance, and gathers the surface operators. Notice the
% line it prints: the inverse is MNE, and the Laplace-Beltrami basis sets the
% finest scale this context can read.

ctx = rheome.flow.context('cfdemo', [], 100);

%%
% The kernel maps the 150 sensors to a three-component current at every
% vertex:

size(ctx.currentKernel)

%% Build the fused kernels
% |rheome.flow.build| turns the kernel into one matrix per flow quantity. Each is
% applied to the sensor data directly; no source time series is formed.

K = rheome.flow.build(ctx);

%% Compute the four maps at one instant
% Take the frame where the sensor signal is strongest, and apply each kernel
% to it. Each map is one matrix-vector product.

[~, t0] = max(vecnorm(ctx.F));
b   = ctx.F(:, t0);
dv  = K.divergence.vertexOperator * b;     % sources +, sinks -
rot = K.curl.vertexOperator * b;           % counter-clockwise +, seen from outside
Phi = K.potential.vertexOperator * b;      % scalar potential
Psi = K.stream.vertexOperator * b;         % stream function

%% Show them
% Draw each map on the cortex, seen from the left.

tl = tiledlayout(2, 2, 'TileSpacing', 'compact');
names = {'Divergence', 'Rotation', 'Phi (scalar potential)', 'Psi (stream function)'};
maps  = {dv, rot, Phi, Psi};
for k = 1:4
    rheome.show.surface(ctx.S, maps{k}, 'Parent', nexttile(tl), 'Title', names{k});
end

%%
% The planted field is a pure rotation, so the rotation and the stream
% function carry it and the divergence is weak. The stream function has one
% extremum, at the vortex.

%% Check the result
% The extremum of Psi should sit at the planted vortex. Measure the distance:

[~, iPeak] = max(abs(Psi));
distMM = 1e3 * norm(ctx.S.Vertices(iPeak, :) - demo.vortexCentre)
assert(distMM < 20)

%% What you have done
% You built the minimum-norm current from a leadfield and a noise
% covariance, fused it with the surface operators into one kernel per
% quantity, and read divergence, rotation, Phi and Psi from a single frame
% of sensor data. Next:
%
% * <howto_flow_maps.html How to compute flow maps> over many frames, in
%   bands, and at a declared scale.
% * <about_resolution_floor.html The resolution floor>: why per-vertex maps
%   like these are not read as physiology, and at which scales the maps
%   can be read.
%
% _Written for Rheome @COMMIT@._

