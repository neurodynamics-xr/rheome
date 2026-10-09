% SCALE  The per-participant driver (rheome.scale.run), its measures and the group reduction.
%
% Functions:
%   rheome.scale.analyses               - The analyses the scale-up covers, and which of them the driver runs yet.
%   rheome.scale.bandsnr                - Per-octave SNR of the rest recording against the subject's own noise run.
%   rheome.scale.cleanspan              - Where a recording is stationary enough to tile: artefact blocks and trimmed edges.
%   rheome.scale.cleanwelch             - Welch PSD of a span, averaging only the segments that touch no artefact block.
%   rheome.scale.coefficients           - Graph-wavelet coefficient envelopes rolled up to cortical tiles: Prognome's MEG input.
%   rheome.scale.importsubject          - One Brainstorm-protocol subject -> the data cache the analyses load from.
%   rheome.scale.locate                 - Find one subject's files inside an unzipped Brainstorm protocol.
%   rheome.scale.measure_bandperiodic   - Table 5 for one subject: per band, rhythm against background and instrument.
%   rheome.scale.measure_bandresolution - Resolution at each octave's OWN measured SNR, for one subject.
%   rheome.scale.measure_catalogue     - The catalogue's patterns through this participant's own head, read by three estimators (MS1 G6, G16).
%   rheome.scale.measure_fusion         - Fused flow kernels against reconstruct-then-differentiate, on one subject's data.
%   rheome.scale.measure_atlas          - MS1 G12 per participant: DK correspondence, atlas events, shared-gauge tensor, connectome roll-up.
%   rheome.scale.measure_geometry       - The atlas geometry checks on one subject's cortex: atom-tile overlap, gauge, roll-up.
%   rheome.scale.measure_grouptrack     - Viterbi tile tracking of the alpha envelope, real against the mode-shift surrogate.
%   rheome.scale.measure_helmholtzbands - Helmholtz-band recovery on the subject's own cortex (geometry only; no MEG).
%   rheome.scale.measure_ownregion      - Readout rule 5 for one subject: which Desikan-Killiany regions own their flux.
%   rheome.scale.measure_inject         - Known movers injected into the own recording, tracked back (MS1 G7, Fig. 8).
%   rheome.scale.measure_trackfactorial - The 0.52 diagnosis as a 2^5 factorial on one participant (MS1 G8, Fig. 9).
%   rheome.scale.refhead                - A cached head without its recording: the factorial's reference head (subject01).
%   rheome.scale.measure_periodicflow   - Apparent alpha flow on the TOTAL and the PERIODIC envelope, over many tiles.
%   rheome.scale.measure_resolution     - The resolution report's headline numbers for one subject.
%   rheome.scale.reduce                 - Per-subject tables -> group distributions, the anchor's place, cohort contrasts.
%   rheome.scale.rows                   - Long-format metric rows: one row per (analysis, metric, band).
%   rheome.scale.run                    - One subject, every ported analysis, compact tables out. The per-subject driver.
%   rheome.scale.sensors                - The channel selection, leadfield, noise covariance and kernel every measure shares.
%
% Author: Diellor Basha, 2026
