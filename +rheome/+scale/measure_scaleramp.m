function [T, X] = measure_scaleramp(name, S, opts)
% SCALE.MEASURE_SCALERAMP  Tracking across scale: a peak whose size ramps 10 -> 40 mm, read by rheome.detect.ridges.
%
%   [T, X] = rheome.scale.measure_scaleramp(name, S)
%   [T, X] = rheome.scale.measure_scaleramp(name, S, SigmaMM=[10 40], NFrames=49, Surfaces=["sphere" "cortex"])
%
% MS1 2.4 (VR1). Can the instrument follow one pattern as it changes SIZE, and read that size in mm? A heat-kernel
% peak sits at one vertex while its width sigma ramps linearly from SigmaMM(1) to SigmaMM(2) over NFrames;
% rheome.detect.ridges reads the field (no noise, no forward model: this calibrates the detector's scale axis,
% not the inverse). Two surfaces:
%   sphere  an ico5 sphere (10242 vertices) with the area of the participant's left hemisphere, its own LBO
%           basis (KSphere modes), built inside the job
%   cortex  the participant's left hemisphere and its cached LBO basis (S.B.L.lbo), the peak at the vertex
%           nearest the hemisphere's area centroid
%
% THE FIELD. f_t = Phi diag(exp(-lambda sigma_t^2 / 2)) Phi' M delta_c, normalised to peak 1 per frame: the
% heat kernel at diffusion time sigma_t^2/2, a Gaussian of width sigma_t (geodesic, m) on a flat patch.
%
% ⭐ THE CALIBRATION FOR A SCALAR PEAK IS gamma = 1, NOT ridges' 1/sqrt(2). The mexhat member
% t lambda e^{-t lambda} meets the peak's spectrum e^{-lambda sigma^2/2} at its centre; with Weyl's
% uniform mode density the response is int t lam e^{-t lam} e^{-lam sigma^2/2} dlam = t / (t + sigma^2/2)^2,
% maximal at t = sigma^2/2, i.e. at the member sigma_m = sqrt(2t) = sigma. ridges' .sigmaTrue
% (= sqrt(2) sigma_m) is the VORTEX calibration (spectrum lambda e^{-lambda sigma^2/2}); it does not apply
% here. The read size reported is the member label sigma_m, refined by a parabola in log sigma through the
% response of the member and its two neighbours at the ridge vertex; the slope against the planted sigma
% should be 1.
%
% ROWS (band = surface): slope and intercept (mm) of read vs planted sigma, ratio_median (read / planted),
% log_err_max (max |log(read / planted)|), continuity (frames on the longest ridge / NFrames), n_ridges,
% n_scales_visited (members the longest ridge visits), centre_err_mm_max (ridge vertex to the planted
% vertex, Euclidean), and the bank's finest usable member. X: one row per surface x frame.
%
% See also: rheome.detect.ridges, rheome.filters.mexhat, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.SigmaMM (1,2) double = [10 40]
        opts.NFrames (1,1) double {mustBeInteger, mustBePositive} = 49   % 0.625 mm a frame: <= 0.7 member at 8 per octave, inside ridges' one-member link
        opts.BankMM (1,2) double = [5 80]           % the member range, wider than the ramp so its ends are interior
        opts.PerOctave (1,1) double = 8
        opts.KSphere (1,1) double = 1000
        opts.SphereSub (1,1) double = 5
        opts.Surfaces string = ["sphere" "cortex"]
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    Hh = S.B.L;  Sc = Hh.S;  Sc.nV = size(Sc.Vertices, 1);
    area = full(sum(Hh.lbo.Mass, 'all'));
    sig = linspace(opts.SigmaMM(1), opts.SigmaMM(2), opts.NFrames) * 1e-3;
    no = log2(opts.BankMM(2) / opts.BankMM(1));
    bank = opts.BankMM(1) * 2.^(0:1/opts.PerOctave:no) * 1e-3;
    T = rheome.scale.rows("scaleramp", strings(0,1), [], "");  X = table();
    for s = opts.Surfaces(:)'
        switch s
            case "sphere"
                [V, F] = rheome.geom.icosphere(opts.SphereSub);  V = V * sqrt(area / (4*pi));
                Sx = struct('Vertices', V, 'Faces', F, 'nV', size(V, 1));
                [L, M] = rheome.operators.laplace_beltrami(V, F);  b = rheome.eigen.modes(L, M, opts.KSphere);
                [~, c] = max(V(:, 3));                                  % the +z pole
            case "cortex"
                Sx = Sc;  b = Hh.lbo;  a = full(sum(b.Mass, 2));
                [~, c] = min(vecnorm(Sx.Vertices - sum(a .* Sx.Vertices) / sum(a), 2, 2));
        end
        r = i_ramp(Sx, b, c, sig, bank);
        T = [T; rheome.scale.rows("scaleramp", ...
             ["slope" "intercept_mm" "ratio_median" "log_err_max" "continuity" "n_ridges" "n_scales_visited" ...
              "centre_err_mm_max" "finest_usable_member_mm" "lambda_max"], ...
             [r.p(1) 1e3*r.p(2) median(r.read ./ sig', 'omitnan') max(abs(log(r.read ./ sig'))) r.cont r.n r.nsc ...
              1e3*max(r.cerr) 1e3*r.fine max(b.Lambda)], ...
             ["ratio" "mm" "ratio" "log" "fraction" "ridges" "members" "mm" "mm" "1/m^2"], s)]; %#ok<AGROW>
        X = [X; table(repmat(s, numel(sig), 1), (1:numel(sig))', 1e3*sig(:), 1e3*r.read, 1e3*r.member, 1e3*r.cerr, r.on, ...
             'VariableNames', {'surface','frame','sigma_planted_mm','sigma_read_mm','sigma_member_mm','centre_err_mm','on_ridge'})]; %#ok<AGROW>
        fprintf('[scaleramp %s] %s: slope %.3f, ratio %.3f, continuity %.2f, %d ridges\n', name, s, r.p(1), ...
                median(r.read ./ sig', 'omitnan'), r.cont, r.n);
    end
end

function r = i_ramp(Sx, b, c, sig, bank)
% plant the ramp, read it by ridges, refine the size on the longest ridge
    d = zeros(Sx.nV, 1);  d(c) = 1;  q = b.Phi' * (b.Mass * d);  lam = b.Lambda(:);
    f = b.Phi * (exp(-lam * sig.^2 / 2) .* q);  f = f ./ max(f, [], 1);
    w = warning('off', 'detect:ridges:truncated');  cw = onCleanup(@() warning(w));
    out = rheome.detect.ridges(f, Sx, b, bank, struct('keepVolume', true));
    sm = out.summary.sigmas;  nT = numel(sig);
    [~, k] = max(arrayfun(@(z) z.nFrames, out.ridges));  R = out.ridges(k);
    r.member = nan(nT, 1);  r.read = nan(nT, 1);  r.cerr = nan(nT, 1);  r.on = false(nT, 1);
    for i = 1:numel(R.t)
        t = R.t(i);  m = R.m(i);  r.on(t) = true;  r.member(t) = sm(m);
        r.read(t) = sm(m);
        if m > 1 && m < numel(sm)                                    % parabola in log sigma through m-1, m, m+1
            y = double(out.W(R.v(i), m-1:m+1, t));  den = y(1) - 2*y(2) + y(3);
            if den < 0, r.read(t) = exp(log(sm(m)) + 0.5 * (y(1) - y(3)) / den * mean(diff(log(sm(m-1:m+1))))); end
        end
        r.cerr(t) = norm(Sx.Vertices(R.v(i), :) - Sx.Vertices(c, :));
    end
    ok = r.on;  r.p = polyfit(sig(ok)', r.read(ok), 1);
    r.cont = R.nFrames / nT;  r.n = out.summary.nRidges;  r.nsc = R.nScales;
    r.fine = min(sm(out.summary.usable));
end

% Author: Diellor Basha, 2026
