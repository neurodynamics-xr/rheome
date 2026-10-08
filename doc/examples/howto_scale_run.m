%% How to run the per-participant measures
% *How-to guide.* |rheome.scale.run| runs the toolbox's per-participant measures on
% one cached participant and writes compact tables; |rheome.scale.reduce| stacks
% many participants' tables into group distributions. This is the route the
% methods paper used for its cohort measures, one participant per task.
% This guide runs it on the synthetic demo study.

name = 'cfdemo';
if ~rheome.load.has(name)                     % the demo study, imported as in the ingest guide
    demo = cfdemo_study();
    rheome.import.surface(name, demo.cortexFile);
    rheome.import.study(name, demo.studyDir, demo.dataName);
    rheome.import.bases(name, 100, 100, struct('overwrite', true));
end
if ~rheome.load.has('cfdemo_noise')           % its noise run
    demo = cfdemo_study();
    rheome.import.surface('cfdemo_noise', demo.cortexFile);
    rheome.import.study('cfdemo_noise', demo.studyDir, demo.noiseName);
end

%% See which measures exist

rheome.scale.analyses()

%%
% |resolution|, |bandsnr|, |bandresolution|, |periodicflow| and |grouptrack|
% run; the others are listed but not yet ported.
%
%% Run one participant
% Name the measures, the noise run, and where the tables go. Without
% |OutDir| they go to |rheome.load.outroot()/scale/<name>|.

outDir = fullfile(rheome.load.outroot(), 'scale');
R = rheome.scale.run(name, Analyses=["resolution" "bandsnr" "bandresolution"], ...
              NoiseStudy="cfdemo_noise", OutDir=outDir);

%%
% A measure that fails is recorded with its error and the others still run,
% so check the status of each:

R.timing(:, {'analysis', 'status', 'seconds'})
assert(all(R.timing.status == "ok"))

%%
% The metrics are one long table: one row per measure, metric and band.

head(R.metrics(R.metrics.analysis == "resolution", {'metric', 'value', 'unit'}), 4)

%%
% The folder |outDir/<name>| holds |metrics.csv|, |bandsnr.csv|,
% |timing.csv| and |provenance.json| (the commit, the MATLAB release and the
% options of the run).
%
% |periodicflow| (apparent alpha flow of the total and the periodic
% envelope, over |FlowTiles| tiles) and |grouptrack| (tile tracking of the
% alpha envelope against a surrogate) need minutes of clean rest; add them
% to |Analyses| on real recordings.
%
%% Reduce many participants
% Run |rheome.scale.run| once per participant into the same |outDir|, then reduce:

G = rheome.scale.reduce(outDir, fullfile(rheome.load.outroot(), 'scale_group'));

%%
% The group folder holds |metrics_long.csv| (every row of every
% participant), |status.csv| and |distribution.csv| (median and quartiles per
% measure, metric and band).
%
%% See also
% <about_scale.html Scale: local to global, fast to slow>,
% <about_resolution_floor.html The resolution floor>,
% <reference/rheome.scale.run.html rheome.scale.run>.
%
% _Written for Rheome @COMMIT@._
