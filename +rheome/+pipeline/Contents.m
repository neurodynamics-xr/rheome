% PIPELINE  A plan as a value: operators chained by their arrows, run against an environment.
%
% The registries say where a field lives (+domain), what it is (+fieldtype) and which
% operator accepts it (rheome.operators.registry). This package turns that into an algorithm: a
% plan is a list of steps whose arrows were checked as they were appended, an environment
% holds the matrices, and running the plan returns the answer together with the record of
% how it was produced.
%
%   rheome.pipeline.env       - the bound resources (surface, kernel, eigenbasis) and their domains
%   rheome.pipeline.start     - begin at a field type on a domain, with a frame axis
%   rheome.pipeline.then      - append an operator; the arrow is checked against the current state
%   rheome.pipeline.window    - append a selection on the frame axis (the bookkeeping step)
%   rheome.pipeline.compile   - can this plan run here, and as what function
%   rheome.pipeline.run       - execute, validating each intermediate, and return the record
%   rheome.pipeline.describe  - the plan as a table, with no environment and nothing run
%
% A plan is metadata and an environment is matrices, deliberately: the same plan runs on
% another subject, and a plan can be printed, compared and stored without either.
%
% Example -- sensors to vorticity coefficients, the chain the flow work actually uses:
%
%   p = rheome.pipeline.start("sensorScalar", e.domains.sensors, Frames=rheome.domain.index(nT, Name="time"));
%   p = rheome.pipeline.window(p, [262 270], Rate=600);
%   p = rheome.pipeline.then(p, "inverse_mne");      % -> ambientVertexWorld on the cortex
%   p = rheome.pipeline.then(p, "curl");             % -> scalarVertex
%   p = rheome.pipeline.then(p, "lb_forward");       % -> coeffScalar on the mode line
%   [g, rec] = rheome.pipeline.run(p, e, F);
%
% See also: rheome.domain.of, rheome.fieldtype.registry, rheome.operators.registry
%
% Author: Diellor Basha, 2026
