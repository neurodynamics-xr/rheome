%% How to import a Brainstorm study
% *How-to guide.* This guide imports one participant's Brainstorm study into
% the toolbox's cache, so that |rheome.flow.context|, |rheome.scale.run| and the |rheome.load.*|
% functions can use it. Brainstorm does not need to be running or installed:
% its files are read from disk as plain MAT-files.
%
%% What the study must contain
% From the subject's anatomy folder:
%
% * the cortex surface the leadfield was computed on (|tess_cortex_*.mat|),
%   with the |Structures| atlas that splits it into left and right
%   hemispheres.
%
%%
% From one study (condition) folder:
%
% * the channel file (|channel_*.mat|);
% * an *unconstrained* MEG head model (|headmodel_surf_os_meg*.mat|, overlapping
%   spheres, three orientations per vertex). A constrained head model is
%   refused: the flow kernels need the full three-component current;
% * the noise covariance (|noisecov_full.mat|);
% * the recording (a |data_*.mat| file, imported, not a raw link).
%
%%
% If the recording was imported into a different condition from the one
% holding the channel file and the head model, pass its absolute path as the
% recording name.
%
%% Choose where the cache goes
% The cache is |rheome.load.root()|. For an installed toolbox it is under your
% |userpath|; in a clone of the repository it is |+data| beside the code.
% To put it elsewhere (a larger disk, node-local scratch), set
% |RHEOME_DATA| before importing:
%
%   setenv('RHEOME_DATA', '/path/to/cache')
%
% This guide uses the synthetic study that comes with the documentation.
% Replace |demo.cortexFile|, |demo.studyDir| and |demo.dataName| with your own
% files.

demo = cfdemo_study();
name = 'cfdemo';                       % the name you will load the participant by

%% Import the cortex, the study and the eigenbases
% Three calls, in this order:

rheome.import.surface(name, demo.cortexFile);             % surface + Laplace-Beltrami operators
rheome.import.study(name, demo.studyDir, demo.dataName);  % channels, leadfield, recording, noise covariance
rheome.import.bases(name, 100, 100, struct('overwrite', true));   % 100 eigenmodes per hemisphere

%%
% Choose the number of eigenmodes per hemisphere (the second argument of
% |rheome.import.bases|) for your mesh. It sets the finest scale every later
% analysis can read, so pin it: the flow analyses of the methods paper use
% 400 per hemisphere for the fused kernels and 1000 for resolution and
% Helmholtz bands, on 10,242-vertex hemispheres. If the warning about
% vertices per half-wavelength appears, the top of the basis is unreliable:
% choose fewer modes or a finer mesh. Several bases can coexist in the
% cache; each is stored under its own K.
%
% |rheome.import.dataset| does the same three steps and also computes the
% relative-Dirac eigenbasis, which only the alternative |'dirac'| source
% estimate needs.
%
%% Import a noise recording (optional)
% The per-octave signal-to-noise analysis of |rheome.scale.run| compares the
% recording with a noise run (an empty-room recording) from the same
% channels. Import it under its own name; only the study part is needed:

rheome.import.surface('cfdemo_noise', demo.cortexFile);
rheome.import.study('cfdemo_noise', demo.studyDir, demo.noiseName);

%% Check what was cached

rheome.load.list()

%%
% |rheome.load.has(name)| tells a script whether a participant is already cached,
% so the import can be skipped on the next run:

assert(rheome.load.has(name))
st = rheome.load.study(name);
assert(size(st.hm.Gain, 2) == 3 * st.hm.nV)    % unconstrained: three columns per vertex

%% Use your own minimum-norm kernel instead (optional)
% |rheome.flow.context| computes its own whitened minimum-norm kernel from the
% leadfield and the noise covariance, so you do not need a Brainstorm
% inverse. To keep the kernel Brainstorm computed, for comparison,
% |rheome.import.inverse(name, studyDir)| caches the study's unconstrained
% kernel-only result and |rheome.load.inverse(name)| reads it back.
%
%% See also
% <howto_flow_maps.html How to compute flow maps>,
% <howto_scale_run.html How to run the per-participant measures>,
% <reference/index.html Function reference> (|import|, |load|).
%
% _Written for Rheome @COMMIT@._
