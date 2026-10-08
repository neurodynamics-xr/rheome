%% How to compute flow maps
% *How-to guide.* This guide computes divergence, rotation and the
% Helmholtz potentials for a cached participant: over a block of frames, as
% envelope and phase in a frequency band, as total quantities, and in
% Helmholtz bands at a declared scale. It assumes the participant is
% imported (see <howto_ingest_brainstorm.html How to import a Brainstorm
% study>); here it is the synthetic demo study.

name = 'cfdemo';
if ~rheome.load.has(name)                     % the demo study, imported as in the ingest guide
    demo = cfdemo_study();
    rheome.import.surface(name, demo.cortexFile);
    rheome.import.study(name, demo.studyDir, demo.dataName);
    rheome.import.bases(name, 100, 100, struct('overwrite', true));
end

%% Build the context and the kernels once per participant
% Pin the number of Laplace-Beltrami modes per hemisphere (the third
% argument): it sets the scale axis, and leaving it empty takes whatever
% basis is cached.

ctx = rheome.flow.context(name, [], 100);
K   = rheome.flow.build(ctx);

%%
% The source estimate is a plain whitened minimum norm (|ctx.Method| is
% |'mne'|). The relative-Dirac estimate is an alternative,
% |rheome.flow.context(name, [], Klbo, struct('Method', 'dirac'))|, after
% |rheome.import.dataset|; it band-limits the current before the inverse.
%
% Each kernel comes in two forms: |coeffOperator| [Ks x C] gives
% Laplace-Beltrami coefficients, |vertexOperator| [V x C] gives a map on
% the vertices.

disp(K.report)

%% Apply a kernel to a block of frames
% A map per frame is a matrix-vector product, so a block of frames is one
% matrix-matrix product:

frames = 1:300;                                    % the first second
dv  = K.divergence.vertexOperator * ctx.F(:, frames);   % [V x 300]
rot = K.curl.vertexOperator       * ctx.F(:, frames);
assert(isequal(size(rot), [ctx.V numel(frames)]))

%%
% For the coefficient form, synthesise the vertex map from the modes:

c   = K.curl.coeffOperator * ctx.F(:, frames);     % [Ks x 300]
rotLB = K.curl.scalarModes * c;                    % the band-limited map

%% Get envelope and phase in a band
% The kernels are real and fixed, so filter and Hilbert-transform the
% sensors, not the sources: the result is the same, at the cost of C
% channels instead of 3V source rows.

[bb, aa] = butter(3, [8 12] / (ctx.sfreq/2));
z   = hilbert(filtfilt(bb, aa, ctx.F.')).';        % analytic sensor signal [C x nT]
Psi = K.stream.vertexOperator * z(:, frames);      % complex stream function
PsiEnvelope = abs(Psi);
PsiPhase    = angle(Psi);

%% Read total quantities
% Energy, enstrophy and helicity are quadratic forms of the sensor vector:

E  = sum(ctx.F(:, frames) .* (K.energy.gram    * ctx.F(:, frames)));   % [1 x 300]
Om = sum(ctx.F(:, frames) .* (K.enstrophy.gram * ctx.F(:, frames)));

%% Read sources and vortices at a declared scale
% Per-vertex divergence and rotation carry the point spread of the inverse
% and are not read as physiology. Read sources and vortices from Helmholtz
% band maps instead: the potentials filtered by a tight Laplace-Beltrami
% wavelet bank. Work on one hemisphere, with its own basis.

B  = rheome.load.bases(name, [], [], 100);
H  = B.L;                                          % left hemisphere
rows = reshape((H.gv(:)' - 1)*3 + (1:3)', [], 1);  % its rows of the current
[~, t0] = max(vecnorm(ctx.F));
J  = ctx.currentKernel(rows, :) * ctx.F(:, t0);    % [3nVh x 1]
hb = rheome.differential.helmholtzbands(J, H.S, H.lbo, Maps=true);
table(hb.wavelengthMM(:), hb.Eirr, hb.Esol, 'VariableNames', {'wavelength_mm', 'source_energy', 'vortex_energy'})

%%
% Choose the band before you look at the data when you read a source;
% for a vortex, the band of largest energy is a valid choice. Then the
% extremum of that band's stream function (|hb.PsiBand(:, m)|) locates the
% vortex. Report the harmonic residual (|hb.harmFrac|) with every
% decomposition.

[~, m] = max(hb.Esol);
[~, iv] = max(abs(hb.PsiBand(:, m)));
fprintf('vortex band %.0f mm, at local vertex %d\n', hb.wavelengthMM(m), iv);
fprintf('harmonic residual %.1f%% of the energy\n', 100*hb.harmFrac);

%%
% Only bands at and above the scales the instrument resolves can be read;
% see <about_resolution_floor.html The resolution floor>.
%
%% See also
% <about_helmholtz_hodge.html The Helmholtz-Hodge split>,
% <howto_figures.html How to make figures>,
% <reference/rheome.flow.build.html rheome.flow.build>,
% <reference/rheome.differential.helmholtzbands.html rheome.differential.helmholtzbands>.
%
% _Written for Rheome @COMMIT@._
