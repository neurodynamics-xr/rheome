function out = simulate(name, varargin)
% FORWARD.SIMULATE  Seed a Dirac current field at a cortical vertex, forward it through the head
% model to sensor timeseries, add a REAL empty-room background, and find the detection threshold.
%
%   out = rheome.forward.simulate('sub01')
%   out = rheome.forward.simulate(name, Family="gabor", Band=[8 16], WavelengthMM=140)
%   out = rheome.forward.simulate(name, Family="resonator", MomentNAm=logspace(-1,2,13))
%
% ⭐ ALMOST NOTHING HERE IS NEW MACHINERY. `rheome.filters.impulse(dbasis, vertex, G, domain, nKeep, dir)`
% already seeds an ORIENTED DIPOLE at a vertex in the Dirac basis and returns [3nV x nT]; the eight
% families in +filters supply the analytic kernel on whichever axis they live on --
% `wave`/`dampedwave`/`diffusion`/`kleingordon` on (lambda, t) with domain 'ts', and
% `travwave`/`resonator`/`gabor`/`stmatern` on (lambda, omega) with domain 'js'. This assembles them
% with the leadfield and a real noise floor and reports the one number that was missing: the source
% moment, in nAm, at which the pattern becomes detectable in the sensor timeseries.
%
% ⚠⚠ THE AMPLITUDE NORMALISER IS RSS, NOT THE NET MOMENT, AND THE DIFFERENCE IS NOT COSMETIC. The
% obvious normaliser is the net dipole moment ||sum_v J(v,t)||, and for a ROTOR that sum is nearly
% zero by construction -- the field cancels vectorially -- so dividing by it multiplies the pattern
% by a huge number and the "1 nAm rotor" that comes out is nothing of the kind. `MomentNAm` here is
% the root-sum-square over all vertices and components at the peak time,
% max_t sqrt(sum_v |J(v,t)|^2), which is finite and positive for every pattern including a rotor.
% ⚠ It is therefore NOT comparable to a published "equivalent dipole moment". Measured here, at the
% same 1 nAm RSS: a single-vertex dipole gives 1.25e-15 T array RMS, the unfiltered mode-space seed
% 9.52e-15 T and a 140 mm mexhat pattern 3.11e-15 T -- so an extended pattern at fixed RSS is
% STRONGER at the sensors than a focal dipole of the same RSS, because its net vector moment is
% larger (2.26 nAm against 1.00). A 140 mm annulus is only 7.9 dB below a focal dipole, not the
% near-silence one might expect from a zero-mean pattern.
%
% ⚠ THE EMPTY ROOM SETS THE RATE AND THE CHANNEL SET, not the subject. The background is 300 Hz with
% 297 channels; the subject is 600 Hz with 300. Resampling a noise recording to match a simulation
% changes its spectrum, so the simulation runs at the EMPTY ROOM's rate instead. And the two channel
% sets are intersected BY NAME -- lining them up by index silently pairs different sensors and gives
% a noise field with the right statistics on the wrong geometry.
% ⭐ Measured here: the 270 good MEG channels carry IDENTICAL names in both files, so the intersect
% is a no-op for them and the index route would in fact have worked. It is kept because the
% non-MEG channels do NOT match (subject "SCLK01" against background "SCLK01-177"), so the
% agreement holds only after the MEG restriction and is not a property of the files.
%
% ⚠ Detection is judged against the empty room's OWN distribution over many windows, not against a
% single noise draw or an analytic formula: the criterion is stated in `.criterion` and is that the
% in-band array RMS exceed the noise-only 95th percentile.
%
% ⚠ What this cannot tell you: whether a detected pattern can be LOCALISED or SIZED. Detection in the
% sensor timeseries is a much weaker question than recovery on the cortex -- plant_scale_omega.m
% recovers location to 43-52 mm and scale not at all below 10 dB. Use this for "is it visible", and
% that for "what comes back".
%
% NAME-VALUE
%   Family       "resonator" (default) | "gabor" | "travwave"   -- on (lambda, omega)
%                "oscillator" | "wave" | "dampedwave" | "diffusion"  -- on (lambda, t)
%   Band         [f_lo f_hi] Hz, default [8 16]; f0 = sqrt(f_lo*f_hi), Q = f0/(f_hi-f_lo)
%   WavelengthMM spatial scale for "gabor" (k0 = 2*pi/lambda), default 140
%   Vertex       seed vertex (local, left hemisphere); default = the median-depth vertex
%   Dir          [3x1] dipole direction, or "normal" (default) | "tangential"
%   MomentNAm    sweep of source strengths, default logspace(-1, 2, 13)
%   Duration     seconds, default 4
%   NoiseTrials  empty-room windows for the null, default 200
%   Hemi         "L" (default)
%   SpinRate     Hz. With a [3nVh x 2] Pattern the two columns are taken as a QUADRATURE PAIR and the
%                field ROTATES: J(t) = cos(2*pi*f*t)*P(:,1) + sin(2*pi*f*t)*P(:,2). ⭐ Because a spin
%                about the normal is multiplication of the complex tangent field by exp(i*theta),
%                z*exp(i*w*t) = cos(w*t)*z + sin(w*t)*(i*z) -- so a rotating vortex is exactly the
%                vortex and its pi/2-spun partner (a source) in quadrature, and the forward stays
%                RANK TWO. Default: the band centre, so the vortex turns once per alpha cycle.
%   Pattern      a prebuilt spatial field [3nVh x 1] (e.g. rheome.flow.seedrotor's .J) used INSTEAD of the
%                mode-space seed. ⭐ It goes through the VERTEX leadfield (rheome.forward.leadfield), not
%                Gmode, so nothing is truncated to the 400-mode basis -- which matters for a rotor,
%                whose fine structure the basis does not hold. The temporal factor is unchanged, so a
%                Pattern run is directly comparable to a mode-space run of the same family.
%   Context      a struct with any of .H (rheome.load.bases hemisphere), .d (rheome.load.dirac), .st, .er
%                (rheome.load.study of subject and empty room) so a sweep loads them once. ⚠ 160 calls
%                reloading these took longer than every forward model in the sweep put together.
%
% OUTPUT (struct out)
%   .threshMomentNAm  ⭐ the detection threshold, solved in closed form (not read off the sweep)
%   .rmsPerNAm        in-band array RMS per nAm -- the whole forward model in one number
%   .moment .snrDB .rmsSignal .rmsNoise .rmsNoiseP95 .detected .threshFromGrid .criterion
%   .B [C x nT] the noiseless sensor timeseries at 1 nAm   .C [nModes x nT] the seed coefficients
%   .Gmode [C x nModes] the leadfield in the Dirac basis -- one column per eigenmode's sensor pattern
%   .channels .fs .vertex .depthMM .family .band .lambdaMM
%
% See also: rheome.filters.impulse, rheome.filters.resonator, rheome.filters.gabor, rheome.forward.reconstruct,
%           rheome.forward.observability, plant_scale_omega, rheome.sensors.generate
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Family', "resonator");
    p.addParameter('Band', [8 16], @(x) numel(x)==2 && x(2)>x(1));
    p.addParameter('WavelengthMM', 140, @isscalar);
    p.addParameter('Vertex', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Dir', "normal");
    p.addParameter('MomentNAm', logspace(-1, 2, 13), @isnumeric);
    p.addParameter('Duration', 4, @isscalar);
    p.addParameter('NoiseTrials', 200, @isscalar);
    p.addParameter('Hemi', "L");
    p.addParameter('Pattern', [], @(x) isempty(x) || isnumeric(x));
    p.addParameter('SpinRate', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Context', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Verbose', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;  vb = o.Verbose;
    fam = lower(string(o.Family));
    f0 = sqrt(prod(o.Band));  Qb = f0/diff(o.Band);

    %% geometry and the hemisphere Dirac basis -- reused from Context when a sweep supplies one
    if ~isempty(o.Context)
        ctx = o.Context;
    else
        ctx = struct();
    end
    if isfield(ctx, 'H'), H = ctx.H; else, Bb = rheome.load.bases(name);  H = Bb.(char(o.Hemi)); end
    S = H.S;  gv = double(H.gv(:))';  nVh = size(S.Vertices,1);
    if isfield(ctx, 'd'), d = ctx.d; else, d = rheome.load.dirac(name); end
    cols = find(d.Hemisphere == double(o.Hemi == "R") + 1);   % 1 = L, 2 = R
    if isempty(cols), cols = 1:d.nModes; end
    rows = reshape((gv' - 1)*4 + (1:4), [], 1);
    db = struct('Phi', d.Phi(rows, cols), 'Mass', d.Mass(rows, rows), ...
                'Lambda', d.Lambda(cols), 'nVert', nVh, 'nModes', numel(cols));

    %% the instrument, and the empty room that sets the rate and the channel set
    if isfield(ctx, 'st'), st = ctx.st; else, st = rheome.load.study(name); end
    if isfield(ctx, 'er'), er = ctx.er; else, er = rheome.load.study('emptyroom'); end
    fs = er.rec.sfreq;  nT = round(o.Duration*fs);
    okS = strcmpi(st.chan.Type,'MEG') & (st.rec.ChannelFlag(:)'==1);
    okE = strcmpi(er.chan.Type,'MEG') & (er.rec.ChannelFlag(:)'==1);
    % ⚠ the empty-room study carries only .Name and .Type -- no .Channel struct, so
    % {chan.Channel.Name} works on the subject and errors on the background.
    nmS = string(st.chan.Name(:))';  nmE = string(er.chan.Name(:))';
    [common, iS, iE] = intersect(nmS(okS), nmE(okE), 'stable');   % ⚠ BY NAME, never by index
    fS = find(okS);  fE = find(okE);  iS = fS(iS);  iE = fE(iE);
    Gain = double(st.hm.Gain(iS, :));
    NF = double(er.rec.F(iE, :));
    Loc = cell2mat(arrayfun(@(c) c.Loc(:,1), st.chan.Channel(iS), 'UniformOutput', false));
    if vb, fprintf('%d channels common to subject and empty room, %g Hz, %.0f s of background\n', ...
            numel(common), fs, size(NF,2)/fs); end

    %% the seed
    dmin = 1e3*min(pdist2(S.Vertices, Loc'), [], 2);
    v = o.Vertex;
    if isempty(v), [~, v] = min(abs(dmin - median(dmin))); end
    if isstring(o.Dir) || ischar(o.Dir)
        n = S.VertNormals(v,:)';
        if lower(string(o.Dir)) == "normal", dir = n;
        else, e = null(n'); dir = e(:,1); end
    else
        dir = o.Dir(:);
    end
    dir = dir/norm(dir);

    %% ⭐ THE LEADFIELD IN THE DIRAC MODE BASIS -- the forward operator this wants
    % The leadfield answers "what sensor pattern does a source delta produce". Composed with the
    % Dirac basis it answers "what sensor pattern does an EIGENMODE produce": Gmode = Gain * Phi,
    % [C x nModes], one column per eigenmode. Seeding is then a choice of coefficients C(t) and the
    % sensors are ONE GEMM, B = Gmode * C. Nothing of size [3nV x nT] is ever built.
    %
    % ⚠⚠ rheome.forward.dirac IS NOT THIS OPERATOR, despite being the obvious candidate and running without
    % complaint. It is the ANALYSIS-convention gain -- what rheome.inverse.dirac consumes -- and it carries
    % a mass weighting: measured against Gain*Phi on the same coefficients the two differ by a factor
    % of 1.23e5, which is 1/(mean vertex area) = 9.1e4 to within the spread of the areas. Used as a
    % synthesis operator it puts the sensor field 100+ dB too low, and nothing in the output looks
    % wrong -- the units are still tesla and the pattern is still the right shape. ⭐ Build the
    % synthesis operator explicitly from the current (imaginary) rows of Phi.
    [Gmode, Psi] = rheome.forward.diracgain(Gain(:, reshape((gv-1)*3 + (1:3)', 1, [])), db);
    if vb, fprintf('leadfield in the Dirac basis: %s; per-mode sensor norm %.2e to %.2e T/coeff\n', ...
            mat2str(size(Gmode)), min(vecnorm(Gmode,2,1)), max(vecnorm(Gmode,2,1))); end

    %% a prebuilt spatial pattern short-circuits the mode-space seed
    usePat = ~isempty(o.Pattern);
    spin2  = false;
    if usePat
        Jp = double(o.Pattern);
        if size(Jp,1) ~= 3*nVh
            Jp = reshape(Jp, 3*nVh, []);
        end
        if size(Jp,2) == 2
            spin2 = true;                                     % a quadrature pair: the field rotates
        elseif size(Jp,2) ~= 1
            error('forward:simulate:pattern', 'Pattern must have 1 or 2 columns, got %d.', size(Jp,2));
        end
        GLv = Gain(:, reshape((gv-1)*3 + (1:3)', 1, []));
        % ⭐ RANK ONE (or TWO), so nothing of size [3nV x nT] is built.
        b1 = GLv * Jp;                                        % [C x 1] or [C x 2]
        rssPat = max(vecnorm(Jp, 2, 1));
    end

    %% the seed, as coefficients on (eigenmode, time) or (eigenmode, frequency)
    lam = db.Lambda(:);
    % ⚠⚠ THE DIRAC SPECTRUM HAS NO LENGTH AXIS, so rheome.filters.mexhat(db.Lambda, (lambda/2pi)^2) is
    % dimensionally meaningless -- and it FAILS QUIETLY: the first sweep returned the same threshold
    % for 267, 189, 140, 95 and 70 mm to three digits, because the filter was landing in the same
    % nearly-flat region of an unrelated spectrum every time. ⭐ Calibrate a wavenumber onto each
    % Dirac mode by projecting its current (imaginary) part onto the LBO and taking the
    % energy-weighted mean wavenumber, then filter on THAT. Same calibration as every other script
    % in this project.
    Lb = H.lbo;
    rxl = (1:nVh)*4 - 2;  ryl = (1:nVh)*4 - 1;  rzl = (1:nVh)*4;
    kb = nan(db.nModes,1);  shr = zeros(db.nModes,1);
    for mm = 1:db.nModes
        e = sum(abs(Lb.Phi'*(Lb.Mass*[db.Phi(rxl,mm) db.Phi(ryl,mm) db.Phi(rzl,mm)])).^2, 2);
        if sum(e) > 0, shr(mm) = sum(e);  kb(mm) = sqrt(sum(e.*Lb.Lambda(:))/sum(e)); end
    end
    wlM = 1e3*2*pi./max(kb, eps);
    okM = isfinite(kb) & kb > 0 & shr > 0.1*max(shr) & wlM < 400;
    k = sqrt(max(lam,0));
    gs = rheome.filters.mexhat(kb.^2, (o.WavelengthMM*1e-3/(2*pi))^2);
    gs(~okM) = 0;                                             % drop the near-null modes (6e9 mm)
    if ~any(gs ~= 0)
        error('forward:simulate:scale', 'no usable Dirac mode near %g mm', o.WavelengthMM);
    end
    % ⚠ the seed's PLACE still has to come from the vertex: a spatial filter alone is a global
    % annulus in the spectrum, identical everywhere. The place enters as the projection of a dipole
    % delta at the vertex onto the modes, which is a row of Phi -- the mode-space image of the seed.
    rv = (v-1)*4 + (2:4);
    a0 = (db.Phi(rv, :).' * dir);                            % [nModes x 1] the seeded dipole
    c0 = gs .* a0;                                           % filtered to the wanted spatial scale
    t  = (0:nT-1)/fs;
    switch fam
        case "resonator"
            % eigenmode-FREQUENCY: a resonance at f0 with bandwidth Q, per mode
            Hf = rheome.filters.resonator(lam, i_omega(fs,nT), f0, Qb);   % [nModes x nF]
            C  = i_fromspectrum(c0 .* Hf, nT);
        case "gabor"
            Hf = rheome.filters.gabor(lam, i_omega(fs,nT), 2*pi/(o.WavelengthMM*1e-3), f0, 1);
            C  = i_fromspectrum(c0 .* Hf, nT);
        case "travwave"
            Hf = rheome.filters.travwave(lam, i_omega(fs,nT), 1.0, 1);
            C  = i_fromspectrum(c0 .* Hf, nT);
        case "oscillator"
            % eigenmode-TIME, the simplest analytic case: every mode oscillates at f0 with a
            % Hann envelope, so the pattern is a standing oscillation of a fixed spatial shape
            C = c0 * (cos(2*pi*f0*t) .* hann(nT)');
        case "wave"
            C = c0 .* cos(k * t);                            % (lambda, t): dispersive wave
        case "dampedwave"
            C = c0 .* (cos(k * t) .* exp(-2*t));
        case "diffusion"
            C = c0 .* exp(-lam * t);
        otherwise, error('forward:simulate:family', 'Unknown family "%s".', fam);
    end
    C = real(C);
    if size(C,2) ~= nT
        error('forward:simulate:length', 'the seed is %d samples, expected %d.', size(C,2), nT);
    end

    %% ⭐ amplitude: normalise in the units the LEADFIELD speaks, A m per vertex
    % ⚠ NOT the coefficient norm. Phi is M-orthonormal, so ||C|| is the MASS-WEIGHTED source norm,
    % which is not sum_v |J_v|^2 and differs from it by the vertex areas. Reconstruct the single
    % peak-time column -- one [4nV x nModes] matvec -- and scale by its per-vertex RSS.
    [~, tp] = max(vecnorm(C, 2, 1));
    if usePat && spin2
        % ⭐ THE ROTATING CASE. Both quadrature components are driven, so the spatial pattern turns
        % rather than merely oscillating: chirality reverses every half turn.
        fs_ = o.SpinRate;  if isempty(fs_), fs_ = f0; end
        tt  = (0:nT-1)/fs;
        env = hann(nT)';
        C   = [cos(2*pi*fs_*tt).*env; sin(2*pi*fs_*tt).*env];   % [2 x nT]
        rssPeak = rssPat * max(vecnorm(C, 2, 1));
    elseif usePat
        % ⚠⚠ THE COLLAPSE TO A SCALAR SERIES MUST KEEP THE SIGN. An earlier version used
        % vecnorm(C,2,1), which is non-negative: for the oscillator family that turns
        % cos(2*pi*f0*t) into |cos(2*pi*f0*t)|, whose fundamental sits at 2*f0 = 22.6 Hz -- OUTSIDE the
        % 8-16 Hz analysis band. The band-pass then removed almost all of it and every Pattern
        % threshold came back inflated by ~1500x. ⭐ The leading left singular vector gives a SIGNED
        % series and is exact for a separable family, where C is rank one.
        [u1, ~, ~] = svds(C, 1);
        tw = real(u1' * C);
        if max(abs(tw)) == 0, error('forward:simulate:null', 'the seed has no time course'); end
        tw = tw / max(abs(tw));
        rssPeak = rssPat * max(abs(tw));
        C = tw;                                               % [1 x nT], signed
    else
        rssPeak = norm(Psi * C(:, tp));                       % per-vertex RSS, in A m
    end
    if ~(rssPeak > 0), error('forward:simulate:null', 'the seeded field is identically zero'); end
    C1 = C / rssPeak * 1e-9;                                 % 1 nAm RSS at the peak time
    if usePat
        B1 = b1 * C1;                                        % [C x 1] x [1 x nT], rank one
    else
        B1 = Gmode * C1;                                     % [C x nT] -- the whole forward model
    end
    % ⭐ verify the mode-space forward against the vertex route on ONE column, since the two differ
    % by a silent mass factor if the wrong operator is used (see the note above).
    % ⚠ A Pattern run IS the vertex route, so there is nothing to cross-check; relErr is 0 by identity
    % rather than by measurement, and the field says which.
    if usePat
        relErr = 0;
    else
        Jchk = Psi * C1(:, tp);
        Bchk = Gain(:, reshape((gv-1)*3 + (1:3)', 1, [])) * Jchk;
        relErr = norm(Bchk - B1(:, tp)) / max(norm(Bchk), realmin);
        if relErr > 1e-10
            error('forward:simulate:route', 'mode-space and vertex forwards disagree by %.2e', relErr);
        end
    end
    J1 = [];                                                 % not materialised; see .C and .Psi
    if vb, fprintf('seed vertex %d, depth %.0f mm, %s, family %s, f0 %.2f Hz, Q %.2f, %.0f mm\n', ...
            v, dmin(v), lower(string(o.Dir)), fam, f0, Qb, o.WavelengthMM); end

    %% the noise-only null, from real empty-room windows
    bp = designfilt('bandpassiir','FilterOrder',8,'HalfPowerFrequency1',o.Band(1), ...
                    'HalfPowerFrequency2',o.Band(2),'SampleRate',fs);
    rng(3);
    starts = randi(size(NF,2)-nT, o.NoiseTrials, 1);
    rmsN = zeros(o.NoiseTrials,1);  W = cell(o.NoiseTrials,1);
    for k = 1:o.NoiseTrials
        w = NF(:, starts(k) + (0:nT-1));
        wf = filtfilt(bp, w.').';
        rmsN(k) = sqrt(mean(wf(:).^2));  W{k} = w;
    end
    thr = prctile(rmsN, 95);
    Bf = filtfilt(bp, B1.').';  rmsS1 = sqrt(mean(Bf(:).^2));

    %% the sweep
    m = o.MomentNAm(:);  nM = numel(m);
    rmsT = zeros(nM,1);  det = false(nM,1);  snr = zeros(nM,1);
    for i = 1:nM
        r = zeros(min(40,o.NoiseTrials),1);
        for k = 1:numel(r)
            y = m(i)*B1 + W{k};
            yf = filtfilt(bp, y.').';
            r(k) = sqrt(mean(yf(:).^2));
        end
        rmsT(i) = median(r);
        det(i)  = rmsT(i) > thr;
        snr(i)  = 20*log10(m(i)*rmsS1/median(rmsN));
    end
    % ⚠⚠ DO NOT READ THE THRESHOLD OFF THE SWEPT GRID. logspace(-1, 2.5, 15) steps by 1.78x, so every
    % true threshold between 1.78 and 3.16 reports as 3.16 -- which is how a sweep over 4 bands and 5
    % wavelengths came back as 3.2 nAm in all 20 cells, a table of the grid rather than of the
    % instrument. ⭐ The threshold has a closed form: the signal RMS is linear in the moment and adds
    % incoherently to the background, so p95^2 = (m*rmsS1)^2 + noise^2 solves exactly.
    k1 = find(det, 1);
    if isempty(k1), thGrid = NaN; else, thGrid = m(k1); end
    ex = thr^2 - median(rmsN)^2;
    if ex > 0, th = sqrt(ex)/rmsS1; else, th = NaN; end

    out = struct('moment', m, 'snrDB', snr, 'rmsSignal', m*rmsS1, 'rmsNoise', median(rmsN), ...
        'rmsNoiseP95', thr, 'rmsTotal', rmsT, 'detected', det, 'threshMomentNAm', th, ...
        'threshFromGrid', thGrid, 'rmsPerNAm', rmsS1, ...
        'criterion', "in-band array RMS above the empty-room 95th percentile over " + ...
                     string(o.NoiseTrials) + " windows", ...
        'B', B1, 'C', C1, 'Gmode', Gmode, 'routeRelErr', relErr, 'J', J1, 'channels', {cellstr(common)}, 'fs', fs, 'vertex', v, ...
        'depthMM', dmin(v), 'family', fam, 'band', o.Band, 'lambdaMM', o.WavelengthMM, ...
        'dir', dir, 'nT', nT, 'route', i_route(usePat), 'spinning', spin2, ...
        'spinRateHz', i_spinrate(spin2, o.SpinRate, f0));
    if vb
        fprintf('1 nAm gives in-band array RMS %.3e T against a noise floor of %.3e T (p95 %.3e)\n', ...
            rmsS1, median(rmsN), thr);
        fprintf('⭐ detection threshold: %.3f nAm (analytic); %.2f from the swept grid\n', th, thGrid);
    end
end

function r = i_spinrate(spin2, sr, f0)
    if ~spin2, r = NaN; elseif isempty(sr), r = f0; else, r = sr; end
end

function r = i_route(usePat)
% ⚠ MATLAB has no string multiplication: string(tf)*"a" + string(~tf)*"b" errors. Third time in this
% session, so it gets a helper rather than another attempt at remembering.
    if usePat, r = "pattern"; else, r = "modes"; end
end

function w = i_omega(fs, nT)
    f = (0:floor(nT/2)) * fs/nT;      % the positive half, as the 'js' domain expects
    w = 2*pi*f;
end

function C = i_fromspectrum(Chat, nT)
% A one-sided spectrum [nModes x nF] to a real time series [nModes x nT], with the conjugate
% symmetry put back. ⚠ Doing this by ifft on the half spectrum alone gives a complex series whose
% real part has half the amplitude and the wrong envelope; the mirror is not optional.
    nF = size(Chat, 2);
    full = zeros(size(Chat,1), nT);
    full(:, 1:nF) = Chat;
    if mod(nT,2) == 0, mir = nF-1:-1:2; else, mir = nF:-1:2; end
    full(:, nT-numel(mir)+1:nT) = conj(Chat(:, mir));
    C = real(ifft(full, nT, 2));
end

% Author: Diellor Basha, 2026
