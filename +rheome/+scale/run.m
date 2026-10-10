function R = run(name, opts)
% SCALE.RUN  One subject, every ported analysis, compact tables out. The per-subject driver.
%
%   R = rheome.scale.run('sub02', ProtoDir='/path/to/protocol', Cohort="norm", Dataset="mydata")
%   R = rheome.scale.run('sub01')                       % already cached: skip the import
%   R = rheome.scale.run(name, ..., Analyses=["resolution" "bandresolution"], OutDir=...)
%
% Writes to OutDir/<name>/ (default rheome.load.outroot()/scale/<name>):
%   metrics.csv       long format: subject, cohort, dataset, analysis, metric, band, value, unit
%   bandsnr.csv       the per-octave SNR the band analyses used
%   flowtiles.csv     periodicflow: one row per tile x (total, periodic)
%   grouptrack.csv    grouptrack: one row per hemisphere x depth x frame rate
%   gaugeflow.csv     periodicflow: the alpha flow in the group gauge (FreeSurfer sphere poles), one
%                     row per chart x band x ico3 sphere patch -- the input of rheome.scale.reducegaugeflow
%   fieldsmooth.csv   fieldsmooth: one row per (band, map, frame) -- Dirichlet wavelength, coherence
%   fieldsmooth_maps.mat   fieldsmooth: the strongest frame's raw and band-limited fields (figure data)
%   eventsensors.mat  eventsensors: the tracked event on the sensor data (figure data)
%   plantfloors.csv, movingvortex.csv, composition.csv, sizeruler.csv, vortexscale.csv, rotation.csv,
%   detection.csv, diracangles.csv   the MS1 group plants (P1, nsp cf-plant-floors): the per-plant /
%                     per-placement tables of rheome.scale.measure_plantfloors, measure_movingvortex and
%                     measure_noisefloor; opt-in, not in the "ported" default
%   bandperiodic.csv  bandperiodic: per (half, band) -- Table 5 on the clean span and its two halves
%   helmholtzbands.csv  helmholtzbands: per hemisphere x band x plant -- Fig. 4A-B on this cortex
%   ownregion.csv     ownregion: per Desikan-Killiany region -- own-region fraction of div and curl
%   geometry.csv      geometry: per hemisphere x level x tile -- atom-tile overlap (Fig. 2A)
%   fusion.csv        fusion: per frame -- fused kernels against reconstruct-then-differentiate
%                     (the five MS1 group P2 tables, nsp rheome-ms1-scale; opt-in, not in "ported")
%   patterns.csv, patternnulls.csv   MS1 G5 / G13 (P3, nsp cf-patterns): rheome.scale.measure_patterns --
%                     detector x band x half rows against the phase-randomised surrogates and the empty room,
%                     and the section 6.3 statistics against their nulls; opt-in, not in the "ported" default
%   catalogue.csv, catalognulls.csv, catalogue_strips.mat   MS1 G6, G3 rule 11, G16 (P4, nsp cf-plants):
%                     rheome.scale.measure_catalogue -- the 19 catalogue plants through the participant's own
%                     gain, whitened MNE and sensor noise (rest, empty room), read by the framework, Brainstorm's
%                     bst_opticalflow and a phase-regression estimator, and the false-propagation nulls;
%                     opt-in, not in the "ported" default. bst_opticalflow needs RHEOME_BRAINSTORM (a Brainstorm checkout)
%   inject.csv, trackfactorial.csv   MS1 G7 / G8 (P5, nsp cf-track): per-injection rows of
%                     rheome.scale.measure_inject and the factorial rows of measure_trackfactorial
%                     (RefHead from RHEOME_REFHEAD); opt-in, not in the "ported" default
%   aperiodic.csv     MS1 G10 (P6, nsp cf-aperiodic): rheome.scale.measure_aperiodic -- per cell (arm x plant x
%                     sigma x speed x signal) the on-patch Viterbi rate against the held-out null, the split's
%                     reduction, band-map errors and along-path speed; opt-in, not in the "ported" default
%   slowosc.csv, slowosc_events.csv   MS1 G15, G16 (P7, nsp cf-slowosc): rheome.scale.measure_slowosc -- sleep
%                     slow oscillations per event x hemisphere x condition (original, reversed, surrogate) x
%                     estimator (framework, bst_of, phasereg, sensorlatency): direction vs A->P, speed; and the
%                     detected events. A sleep EEG subject: Modality="EEG", and Edf/Events with ProtoDir (the
%                     template forward model) to import it (rheome.scale.importeeg); opt-in, not "ported"
%   correspondence.csv, atlasevents.csv, gaugetensor.csv, connectome.csv   MS1 G12 (P8, nsp cf-atlas):
%                     rheome.scale.measure_atlas -- DK agreement through the sphere vs raw coordinates and the
%                     frame vs the shared meridian; alpha atlas events at depth <= 3 vs phase-randomised
%                     surrogates; the alpha orientation tensor per depth-3 group tile in the shared gauge;
%                     Destrieux -> DK vs the dyadic roll-up and connectome-wavelet widths over a gamma sweep.
%                     Opt-in. Needs a template cortex (RHEOME_TEMPLATE or RHEOME_BRAINSTORM's ICBM152)
%   multimodal.csv, multimodal_assoc.csv, multimodal_connectome.mat   the MS1 multimodal example (ea515f28):
%                     rheome.scale.measure_multimodal -- MEG band power, PET SUVR and the fibre-degree density
%                     through one wavelet bank on one dyadic tile ladder, every level an exact roll-up; the
%                     descriptive cross-tile rho; the tile connectome per level. Opt-in. PET and fibres come
%                     from the protocol (pet.mat, fibers_subject.mat in the cache)
%   timing.csv        analysis, seconds, peak_rss_GB, status, message
%   rheome_coeffs.mat Analyses="coefficients" only: graph-wavelet envelopes per tile x scale x time,
%                     Prognome's contract (rheome.scale.coefficients)
%   provenance.json   commit (+ dirty flag), MATLAB release, host, options, located inputs, the
%                     kernels the protocol ships, import durations
% No figures. An analysis that errors is recorded as status "error" with its message and the
% others still run, so one bad subject never costs the array its other numbers.
%
% ⭐ On a cluster, set RHEOME_DATA to node-local disk before calling: the import writes the
% subject's cache there (rheome.load.root), and nothing lands in the checkout.
%
% Ported analyses (the reports' names): resolution, bandsnr, bandresolution, periodicflow
% (apparent flow of the total and the periodic alpha envelope over FlowTiles tiles; its "total"
% rows are the flowmap pipeline's speed, curl and div), grouptrack, fieldsmooth (per-vertex vs
% graph-wavelet band-limited current, div and curl), eventsensors (needs grouptrack first: the
% best Viterbi path's samples and channels on the sensor data), and the MS1 group P2 measures:
% bandperiodic (per-band periodic fraction, oscillation SNR, aperiodic fit, floors, IAF, split
% halves), helmholtzbands (planted-band recovery on the cortex), ownregion (readout rule 5), geometry
% (atom-tile overlap, gauge, roll-up exactness) and fusion (fused-kernel exactness). The others --
% flowmap (div/curl maps), inject, vortex, sensorwavelet -- are listed in
% best Viterbi path's samples and channels on the sensor data). The others --
% flowmap (div/curl maps), vortex, sensorwavelet -- are listed in
% rheome.scale.analyses with status "not_ported" until each has a figure-free measure.
%
% See also: rheome.scale.importsubject, rheome.scale.measure_resolution, rheome.scale.reduce, rheome.scale.analyses
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        opts.ProtoDir char = ''
        opts.Sub char = ''
        opts.Cohort string = ""
        opts.Dataset string = ""
        opts.Analyses string = rheome.scale.analyses("ported")
        opts.OutDir char = ''
        opts.DurationS (1,1) double = 300
        opts.NoiseStudy char = ''      % a cached study whose .rec is the noise (e.g. 'emptyroom')
        opts.FlowTiles (1,1) double = 12   % periodicflow: evenly spaced 2 s tiles per subject
        opts.Modality (1,1) string = "MEG" % "EEG": sleep EEG (rheome.scale.importeeg), average reference
        opts.Edf char = ''                 % with Events: import from EEG-BIDS + ProtoDir's forward model
        opts.Events char = ''
    end
    if isempty(opts.OutDir), opts.OutDir = fullfile(rheome.load.outroot(), 'scale'); end
    od = fullfile(opts.OutDir, name);  if ~exist(od, 'dir'), mkdir(od); end
    P = struct('subject', name, 'cohort', opts.Cohort, 'dataset', opts.Dataset, ...
               'started', char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''')), ...
               'commit', i_git('rev-parse HEAD'), 'dirty', ~isempty(i_git('status --porcelain --untracked-files=no')), ...
               'branch', i_git('rev-parse --abbrev-ref HEAD'), 'matlab', version, 'host', i_host(), ...
               'dataRoot', rheome.load.root(), 'options', opts);

    tim = table('Size',[0 5],'VariableTypes',{'string','double','double','string','string'}, ...
                'VariableNames',{'analysis','seconds','peak_rss_GB','status','message'});
    if ~isempty(opts.ProtoDir)
        t0 = tic;
        try
            if isempty(opts.Edf)
                P.import = rheome.scale.importsubject(name, opts.ProtoDir, DurationS=opts.DurationS, Sub=opts.Sub);
            else
                P.import = rheome.scale.importeeg(name, opts.ProtoDir, opts.Edf, opts.Events);
            end
            tim = [tim; {"import", toc(t0), i_rss(), "ok", ""}];
        catch e
            tim = [tim; {"import", toc(t0), i_rss(), "error", string(e.message)}];
            i_write(od, P, table(), tim, table());  R = struct('provenance', P, 'timing', tim);  return
        end
    end

    M = table();  snr = table();  S = [];  X = struct();  Best = [];  ctx = [];  Kf = [];
    for a = opts.Analyses(:)'
        t0 = tic;  msg = "";  st = "ok";
        try
            if isempty(S) && ~ismember(a, ["bandsnr" "helmholtzbands"]), S = rheome.scale.sensors(name, Modality=opts.Modality); end
            if isempty(Kf) && ismember(a, ["ownregion" "fusion"])          % the fused kernels, built once
                ctx = rheome.flow.context(name);  Kf = rheome.flow.build(ctx);
            end
            switch a
                case "resolution",     M = [M; rheome.scale.measure_resolution(name, S)]; %#ok<AGROW>
                case "bandsnr"
                    if isempty(opts.NoiseStudy), snr = rheome.scale.bandsnr(name);
                    else, nz = rheome.load.study(opts.NoiseStudy); snr = rheome.scale.bandsnr(name, i_noiserec(nz)); end
                    M = [M; rheome.scale.rows("bandsnr", repmat("snr_dB",height(snr),1), snr.snr_dB, "dB", snr.band)]; %#ok<AGROW>
                case "bandresolution"
                    if isempty(snr), error('scale:run:order', 'bandresolution needs bandsnr first'); end
                    M = [M; rheome.scale.measure_bandresolution(name, S, snr)]; %#ok<AGROW>
                case "periodicflow"
                    [Tf, Xf, Gf] = rheome.scale.measure_periodicflow(name, S, NumTiles=opts.FlowTiles);
                    M = [M; Tf];  X.flowtiles = Xf; %#ok<AGROW>
                    if ~isempty(Gf)
                        X.gaugeflow = [table(repmat(string(name),height(Gf),1), repmat(opts.Cohort,height(Gf),1), ...
                            repmat(opts.Dataset,height(Gf),1), 'VariableNames', {'subject','cohort','dataset'}), Gf];
                    end
                case "grouptrack"
                    [Tg, Xg, Best] = rheome.scale.measure_grouptrack(name, S);
                    M = [M; Tg];  X.grouptrack = Xg; %#ok<AGROW>
                case "fieldsmooth"
                    [Ts, Xs, maps] = rheome.scale.measure_fieldsmooth(name, S);
                    M = [M; Ts];  X.fieldsmooth = Xs; %#ok<AGROW>
                    save(fullfile(od, 'fieldsmooth_maps.mat'), 'maps', '-v7.3');
                case "eventsensors"
                    if isempty(Best), error('scale:run:order', 'eventsensors needs grouptrack first'); end
                    [Te, ev] = rheome.scale.measure_eventsensors(name, S, Best);
                    M = [M; Te]; %#ok<AGROW>
                    save(fullfile(od, 'eventsensors.mat'), 'ev', '-v7.3');
                case "plantfloors"             % MS1 G1/G9: Helmholtz-band floors through this MEG
                    [Tp, X.plantfloors] = rheome.scale.measure_plantfloors(name, S);  M = [M; Tp]; %#ok<AGROW>
                case "movingvortex"            % MS1 G11: a moving vortex planted and read back
                    [Tv, X.movingvortex] = rheome.scale.measure_movingvortex(name, S);  M = [M; Tv]; %#ok<AGROW>
                case {"composition" "sizeruler" "vortexscale" "rotation" "detection" "diracangles"}   % MS1 G2
                    [Tn, X.(a)] = rheome.scale.measure_noisefloor(name, S, a);  M = [M; Tn]; %#ok<AGROW>
                case "bandperiodic"
                    nz = [];  if ~isempty(opts.NoiseStudy), nz = i_noiserec(rheome.load.study(opts.NoiseStudy)); end
                    [Tb, X.bandperiodic] = rheome.scale.measure_bandperiodic(name, S, Noise=nz);
                    M = [M; Tb]; %#ok<AGROW>
                case "helmholtzbands"
                    [Th, X.helmholtzbands] = rheome.scale.measure_helmholtzbands(name);
                    M = [M; Th]; %#ok<AGROW>
                case "ownregion"
                    [To, X.ownregion] = rheome.scale.measure_ownregion(name, S, Kf);
                    M = [M; To]; %#ok<AGROW>
                case "geometry"
                    [Tq, X.geometry] = rheome.scale.measure_geometry(name, S);
                    M = [M; Tq]; %#ok<AGROW>
                case "fusion"
                    [Tz, X.fusion] = rheome.scale.measure_fusion(name, ctx, Kf);
                    M = [M; Tz]; %#ok<AGROW>
                case {"patterns" "patternnulls"}   % MS1 G5, G13: catalogue detectors and section 6.3 nulls
                    [Tq, X.(a)] = rheome.scale.measure_patterns(name, S, a);  M = [M; Tq]; %#ok<AGROW>
                case {"catalogue" "catalognulls"}   % MS1 G6, G3 rule 11, G16: catalogue plants and comparators
                    [Tc, X.(a), strips] = rheome.scale.measure_catalogue(name, S, a);  M = [M; Tc]; %#ok<AGROW>
                    if a == "catalogue", save(fullfile(od, 'catalogue_strips.mat'), '-struct', 'strips', '-v7.3'); end
                case "inject"                  % MS1 G7 (Fig. 8): movers injected into the own recording, tracked
                    [Ti, X.inject] = rheome.scale.measure_inject(name, S);  M = [M; Ti]; %#ok<AGROW>
                case "trackfactorial"          % MS1 G8 (Fig. 9): the 0.52 diagnosis, 2^5 factorial
                    [Tf, X.trackfactorial] = rheome.scale.measure_trackfactorial(name, S);  M = [M; Tf]; %#ok<AGROW>
                case "aperiodic"               % MS1 G10: a moving 1/f change, threshold and split selectivity
                    [Ta, X.aperiodic] = rheome.scale.measure_aperiodic(name, S);  M = [M; Ta]; %#ok<AGROW>
                case "slowosc"                 % MS1 G15, G16: sleep slow oscillations, the positive control
                    [Tq, X.slowosc, X.slowosc_events] = rheome.scale.measure_slowosc(name, S);  M = [M; Tq]; %#ok<AGROW>
                case {"correspondence" "atlasevents" "gaugetensor" "connectome"}   % MS1 G12 (P8, nsp cf-atlas)
                    [Ta, X.(a)] = rheome.scale.measure_atlas(name, S, a);  M = [M; Ta]; %#ok<AGROW>
                case "multimodal"              % MS1 multimodal example (ea515f28): MEG + PET + fibres on one ladder
                    [Tm, Xm, Cn] = rheome.scale.measure_multimodal(name, S);  M = [M; Tm]; %#ok<AGROW>
                    X.multimodal = Xm.tiles;  X.multimodal_assoc = Xm.assoc;
                    save(fullfile(od, 'multimodal_connectome.mat'), '-struct', 'Cn', '-v7');
                case "coefficients"            % Prognome's MEG input: opt-in, not in the "ported" default
                    st0 = rheome.load.study(name);
                    Cf = rheome.scale.coefficients(S.B, S.Res.ImagingKernel, double(st0.rec.F(S.iSel,:)), st0.rec.sfreq);
                    clear st0
                    save(fullfile(od, 'rheome_coeffs.mat'), '-struct', 'Cf', '-v7');
                    M = [M; rheome.scale.rows("coefficients", ["n_tiles" "n_scales" "n_samples"], ...
                                              size(Cf.W, 1:3), ["tiles" "scales" "samples"])]; %#ok<AGROW>
                    clear Cf
                otherwise, st = "not_ported";
            end
        catch e
            st = "error";  msg = string(e.identifier) + ": " + string(e.message);
        end
        tim = [tim; {a, toc(t0), i_rss(), st, msg}]; %#ok<AGROW>
        fprintf('[rheome.scale.run %s] %-16s %-10s %7.1f s  %s\n', name, a, st, toc(t0), msg);
    end
    if ~isempty(M)
        n = height(M);
        M = [table(repmat(string(name),n,1), repmat(opts.Cohort,n,1), repmat(opts.Dataset,n,1), ...
                   'VariableNames', {'subject','cohort','dataset'}), M];
    end
    P.finished = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
    i_write(od, P, M, tim, snr, X);
    R = struct('provenance', P, 'metrics', M, 'timing', tim, 'bandsnr', snr, 'tables', X);
end

function rec = i_noiserec(nz)
% a cached noise study -> the rec struct rheome.scale.bandsnr expects (channel names from its chan)
    rec = nz.rec;
    if ~isfield(rec, 'ChannelName') || isempty(rec.ChannelName)
        rec.ChannelName = cellstr(nz.chan.Name);          % rheome.io.read.channel and the ER cache both carry .Name
    end
    if ~isfield(rec, 'sfreq') || isempty(rec.sfreq), rec.sfreq = 1/mean(diff(rec.Time)); end
end

function i_write(od, P, M, tim, snr, X)
    if nargin < 6, X = struct(); end
    if ~isempty(M), writetable(M, fullfile(od, 'metrics.csv')); end
    for f = string(fieldnames(X))'                  % flowtiles.csv, grouptrack.csv
        if ~isempty(X.(f)), writetable(X.(f), fullfile(od, f + ".csv")); end
    end
    if ~isempty(snr), writetable(snr, fullfile(od, 'bandsnr.csv')); end
    writetable(tim, fullfile(od, 'timing.csv'));
    fid = fopen(fullfile(od, 'provenance.json'), 'w');
    fprintf(fid, '%s', jsonencode(P, PrettyPrint=true));  fclose(fid);
end

function s = i_git(args)
    here = fileparts(fileparts(mfilename('fullpath')));
    [rc, s] = system(sprintf('git -C "%s" %s', here, args));
    s = strtrim(s);  if rc ~= 0, s = ''; end
end

function h = i_host()
    [~, h] = system('hostname');  h = strtrim(h);
end

function g = i_rss()
% peak resident set of this MATLAB, GB (VmHWM on Linux; ps on macOS gives the current RSS)
    g = NaN;
    if isfile('/proc/self/status')
        t = fileread('/proc/self/status');  k = regexp(t, 'VmHWM:\s*(\d+)', 'tokens', 'once');
        if ~isempty(k), g = str2double(k{1}) / 1048576; end
    else
        [rc, s] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
        if rc == 0, g = str2double(strtrim(s)) / 1048576; end
    end
end

% Author: Diellor Basha, 2026
