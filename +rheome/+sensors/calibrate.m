function c = calibrate(G, varargin)
% SENSORS.CALIBRATE  Give a dimensionless graph spectrum its length units.
%
%   c = rheome.sensors.calibrate(G)
%
% A Gaussian-kernel graph Laplacian has dimensionless lambda. It SELF-CALIBRATES: for a
% translation-invariant sampling,
%
%   lambda(k) = sum_j w_ij (1 - cos(k . d_ij))  ~  alpha k^2   for small k
%
% and averaging (k.d)^2 over directions in `dim` dimensions gives the CLOSED FORM
%
%   alpha = (1 / (2*dim)) * sum_j w_ij |d_ij|^2                        [metres^2]
%
% ⭐ TWO INDEPENDENT ROUTES, AND THEY MUST AGREE. The closed form above, and a REGRESSION:
% plant cos(k.x) over a range of |k| and directions, take the Rayleigh quotient
% (u'Lu)/(u'u), and fit against k^2 over the small-k region. One route alone is a
% definition; two that agree is a measurement, and .agree is the number to report.
%
% ⭐ THE DEPARTURE IS NOT AN ERROR TERM. lambda_hat(k) falls BELOW alpha*k^2 as the
% wavelength approaches 2*pitch. That curve is the graph's own dispersion relation and it
% marks where the array stops measuring wavelength correctly -- the resolution floor,
% DERIVED rather than asserted. It is returned (.kGrid, .lambdaHat) precisely so a demo can
% show it instead of quoting a rule of thumb.
%
% ⚠ THE ALIAS FLOOR IS NOT THE ACCURACY FLOOR, and .lambdaUsable is the difference.
% lambda ~ alpha k^2 is a SMALL-k expansion: it is already 10% low at roughly 5-6 sensors
% per cycle, long before the 2-sensor alias limit at 2*pitch. Between lambdaUsable and
% 2*pitch a wavelength is still distinguishable but no longer accurately recovered by the
% quadratic model. Any wavelength claim rests on lambdaUsable, not on the alias floor.
% .kUsable/.lambdaUsable are NaN only when the departure is SUSTAINED across the two
% longest-wavelength probed points, not merely crossed once at the first -- an isolated
% departure at point 1 alone does not trigger it (see the sustained-crossing rule below).
% The usable band is therefore "no SUSTAINED departure", not "no departure at all": an
% isolated point inside the band can still exceed 10% between two points that hold below it.
%
% ⚠ NO FIT WINDOW MEANS ROUTE 2 IS NOT TRUSTWORTHY, AND .agree SAYS SO ONLY IF YOU CHECK
% .hasFitWindow FIRST. The small-k regression needs >=3 wavenumbers below FitFraction*kMax
% that are ALSO >=1.5 cycles across the aperture (route 2's own long-wavelength margin,
% above). hasFitWindow is exactly nnz(kGrid <= FitFraction*kMax) >= 3 -- a MEASURED count
% on the actual (logarithmic) wavenumber grid, not a closed-form threshold. A naive
% continuous-k argument suggests Aperture/Pitch >= 3/FitFraction (12 with the defaults),
% but kGrid is discrete and log-spaced, so its point spacing matters too: with the default
% NumK = 24 the measured boundary is Aperture/Pitch ~ 13.7 (it moves if NumK does), not 12.
% A 10x10-at-400um array (Aperture/Pitch = 12.7, e.g. rheome.sensors.grid('utah')) clears the
% naive bar of 12 and STILL has no fit window -- only 1 of the needed 3 wavenumbers lands
% inside it. An 8x8-at-10mm array (Aperture/Pitch ~ 9.9) is not even close: no wavenumber
% is both long enough to sit in the fit window and short enough to fit >=1.5 cycles across
% the array, because the array itself is too small to have a small-k regime. calibrate
% DOES NOT invent one: it warns ('sensors:calibrate:narrowFitWindow'), falls back to the 3
% wavenumbers closest to kMin, and reports .hasFitWindow = false and .nFit < the intended
% window size. It does NOT swap .alpha for the closed form -- that would hide which
% estimator produced the number, which is the same class of problem as answering as if a
% fit window existed. On such an array, .agree is not "two routes agreeing"; it is a good
% estimate (the closed form) checked against one already known to be biased (the
% fallback regression) -- a real instrument limitation, not a bug to average away.
% .kFit reports what was ACTUALLY fitted (max of the selected wavenumbers); the nominal
% target is kept separately as .kFitRequested so the two can be compared.
%
% ⭐ SPEED IS A MULTIPLICATION, NOT A DIVISION. rheome.dynamics.dispersion fits omega =
% c_graph*sqrt(lambda); since sqrt(lambda) = sqrt(alpha)*k, the PHYSICAL speed is
% c_phys = c_graph * speedScale, with speedScale = sqrt(alpha). Units check:
% c_graph is [1/s] (a rate against the dimensionless sqrt(lambda)), speedScale is
% [m], and their product is [1/s]*[m] = [m/s]. Dividing by speedScale instead of
% multiplying makes every speed in the project wrong by a factor of alpha.
%
% ⚠ .alphaSpread IS AN ERROR BAR, NOT DIAGNOSTICS. The closed form assumes every sensor sees
% the same neighbourhood. On an irregular array -- a helmet, a cap -- it does not, and the
% per-vertex spread is a genuine uncertainty on every speed derived from that array. Report
% it WITH the speed. At sensor-array sizes this spread is dominated by BOUNDARY FRACTION,
% not by geometric irregularity as such -- a quasi-uniform cap with a small boundary
% fraction can have a smaller alphaSpread than a small regular lattice with a large one.
%
% ⚠ THE NORMALIZED LAPLACIAN IS REFUSED. D^-1/2 W D^-1/2 rescales each vertex by its own
% degree, so sum_j w_ij |d_ij|^2 is no longer the quantity in the small-k expansion and the
% closed form is simply not that operator's constant.
%
% INPUTS:
%   G  a rheome.sensors.graph struct (combinatorial)
%   'NumK'         wavenumbers on the probe grid (default 24)
%   'Directions'   [nD x 3] unit vectors (default: 8 in the array plane, 1 along a chain)
%   'FitFraction'  fit lambda_hat over k <= FitFraction * (pi/pitch) (default 0.25)
% OUTPUT:
%   c  .alpha .alphaClosed .alphaPerVertex .alphaSpread .agree
%      .kGrid .lambdaHat .kFit .kFitRequested .nFit .hasFitWindow
%      .kUsable .lambdaUsable .speedScale .Directions
%
% See also: rheome.sensors.graph, rheome.sensors.window, rheome.sensors.sample
%
% Author: Diellor Basha, 2026

    if ~strcmp(G.Laplacian, 'combinatorial')
        error('sensors:calibrate:laplacian', ...
            ['alpha is defined for the COMBINATORIAL Laplacian. ''%s'' rescales each ' ...
             'vertex by its degree, so sum_j w_ij |d_ij|^2 is not its constant.'], G.Laplacian);
    end

    arr = G.Array;
    dim = arr.Dim;

    % kMin (below) is 1.5 cycles across the aperture and kMax is the pi/Pitch alias floor;
    % kMin > kMax -- logspace running backwards, and the fallback fitting the three WORST
    % points -- happens whenever Aperture < 3*Pitch. Refuse outright rather than return
    % nonsense from an array too small to have a wavenumber grid at all.
    if arr.Aperture < 3 * arr.Pitch
        error('sensors:calibrate:aperture', ...
            ['%s: Aperture/Pitch = %.2f is too small for route 2''s wavenumber grid -- ' ...
             'the long-wavelength margin (kMin, 1.5 cycles across the aperture) would ' ...
             'exceed the alias floor (kMax, 2 samples per cycle). Needs Aperture >= ' ...
             '3*Pitch; got aperture %.3g m at pitch %.3g m.'], ...
            arr.Name, arr.Aperture / arr.Pitch, arr.Aperture, arr.Pitch);
    end

    p = inputParser;
    p.addParameter('NumK',        24);
    p.addParameter('Directions',  []);
    p.addParameter('FitFraction', 0.25);
    p.parse(varargin{:});
    o = p.Results;

    % ---- route 1: the closed form, per vertex ----
    alphaPerVertex = full(sum(G.W .* (G.D.^2), 2)) / (2*dim);
    alphaClosed    = mean(alphaPerVertex);
    alphaSpread    = std(alphaPerVertex) / max(alphaClosed, eps);

    % ---- probe directions ----
    dirs = o.Directions;
    if isempty(dirs)
        dirs = i_directions(arr);
    end
    dirs = dirs ./ vecnorm(dirs, 2, 2);

    % ---- route 2: plant plane waves, take the Rayleigh quotient ----
    % kMin is NOT 2*pi/Aperture. At exactly one wavelength across the aperture the planted
    % wave is a single, nearly-DC cycle over a FINITE, non-periodic domain, and the Rayleigh
    % quotient measured on the actual graph departs from the true small-k value by 20-30%
    % (checked against the exact infinite-lattice sum, which tracks alpha*k^2 to <0.2% at
    % that same k -- so the departure is a finite-domain sampling artifact of route 2 itself,
    % not lattice dispersion). A 1.5-cycle margin -- half the margin the alias floor keeps
    % at the short-wavelength end -- puts kGrid's first point back on the quadratic model to
    % ~1%, on both a fine 32x32 lattice and a coarse 8x8 preset, and leaves more of the
    % small-k grid available for the fit below than a larger margin would.
    kMin = 1.5 * (2*pi / arr.Aperture);
    kMax = pi   / arr.Pitch;                       % the alias floor, in wavenumber
    kGrid = logspace(log10(kMin), log10(kMax), o.NumK);

    lambdaHat = zeros(1, numel(kGrid));
    for i = 1:numel(kGrid)
        acc = 0;
        for d = 1:size(dirs, 1)
            kv = kGrid(i) * dirs(d, :);
            u  = rheome.sensors.sample(arr, @(P) cos(P * kv.'));
            u  = u - mean(u);                      % the constant mode carries no wavenumber
            nu = u.' * u;
            if nu > eps
                acc = acc + (u.' * (G.L * u)) / nu;
            end
        end
        lambdaHat(i) = acc / size(dirs, 1);
    end

    % ---- fit alpha over the small-k region only ----
    kFitRequested = o.FitFraction * kMax;
    selWindow     = kGrid <= kFitRequested;
    hasFitWindow  = nnz(selWindow) >= 3;
    if hasFitWindow
        sel = selWindow;
    else
        % NO SMALL-k REGIME EXISTS ON THIS ARRAY. Falling back to the 3 wavenumbers
        % closest to kMin is the least-bad option, but it is NOT the fit window the
        % caller asked for, and .alpha is NOT swapped for the closed form -- see the
        % header. Warn with exactly what was measured instead.
        sel = 1:min(3, numel(kGrid));
        spc = (2*pi ./ kGrid(sel)) / arr.Pitch;    % sensors per cycle actually fitted
        warning('sensors:calibrate:narrowFitWindow', ...
            ['%s: Aperture/Pitch = %.1f puts only %d of the needed 3 wavenumbers inside ' ...
             'the small-k fit window (k <= %.3g rad/m; the boundary is ~13.7 with the ' ...
             'default NumK = 24, not a fixed multiple of FitFraction, since it depends ' ...
             'on NumK). Falling back to the %d wavenumbers closest to kMin, at %.1f-%.1f ' ...
             'sensors/cycle -- already past the dispersive floor on this array. ' ...
             '.alpha remains route 2''s regression value, not the closed form.'], ...
            arr.Name, arr.Aperture / arr.Pitch, nnz(selWindow), kFitRequested, numel(sel), ...
            min(spc), max(spc));
    end
    nFit     = numel(kGrid(sel));
    kFit     = max(kGrid(sel));                    % what was ACTUALLY fitted
    kk       = kGrid(sel).^2;
    alphaEmp = (kk * lambdaHat(sel).') / (kk * kk.');    % least squares through the origin

    % ---- where the quadratic model stops being the truth ----
    % ⚠ THE FIRST CROSSING IS NOT THE FLOOR ON A 1D ARRAY. lambda_hat(k) is averaged over
    % the probe directions and a chain admits exactly ONE (see i_directions), so its curve
    % keeps route 2's finite-domain artifact undamped and an ISOLATED point can cross 10%
    % many sensors-per-cycle before the departure is sustained. Require the crossing to
    % HOLD at the next grid point too: the floor is where the model stops being the truth,
    % not where it flickers.
    model = alphaEmp * kGrid.^2;
    over  = abs(lambdaHat ./ max(model, eps) - 1) > 0.10;
    bad   = find(over(1:end-1) & over(2:end), 1, 'first');
    if isempty(bad)
        kUsable = kMax;
    elseif bad == 1
        % The model has already failed at the LONGEST wavelength probed -- there is no
        % usable band on this array, and returning kGrid(1) would read as a resolution
        % claim it isn't.
        kUsable = NaN;
    else
        kUsable = kGrid(bad - 1);
    end

    c = struct('alpha', alphaEmp, 'alphaClosed', alphaClosed, ...
               'alphaPerVertex', alphaPerVertex, 'alphaSpread', alphaSpread, ...
               'agree', abs(alphaEmp - alphaClosed) / max(alphaClosed, eps), ...
               'kGrid', kGrid, 'lambdaHat', lambdaHat, ...
               'kFit', kFit, 'kFitRequested', kFitRequested, ...
               'nFit', nFit, 'hasFitWindow', hasFitWindow, ...
               'kUsable', kUsable, 'lambdaUsable', 2*pi/kUsable, ...
               'speedScale', sqrt(alphaEmp), 'Directions', dirs);
end

function dirs = i_directions(arr)
% Probe directions IN THE ARRAY'S OWN PLANE. A chain admits exactly one, which is the same
% fact as its inability to resolve direction -- stated by the geometry, not special-cased.
    Pc = arr.Pos - mean(arr.Pos, 1);
    [~, ~, V] = svd(Pc, 'econ');
    if arr.Dim < 2
        dirs = V(:,1).';
        return;
    end
    th = (0:7).' * (pi/8);                       % 8 directions over a half-turn
    dirs = cos(th) * V(:,1).' + sin(th) * V(:,2).';
end

% Author: Diellor Basha, 2026
