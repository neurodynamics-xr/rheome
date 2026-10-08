function out = blobscale(x, basis, opts)
% DETECT.BLOBSCALE  The characteristic SIZE of structure in a cortical field, by scale selection.
%
%   B = rheome.load.bases('sub01');  basis = B.L.lbo;        % cached: 1000 modes, no solve
%   out = rheome.detect.blobscale(x, basis)
%   out = rheome.detect.blobscale(x, basis, Kernel="mexhat", Wavelengths=[40 500]*1e-3)
%   out = rheome.detect.blobscale(curl, basis, MinWavelength=0.122)     % floor at the array's resolution
%
% ⚠⚠⚠ .wavelength IS A MAP [nV x nT], NOT THE SIZE OF THE BLOB. It is the scale selected AT
% EACH VERTEX, so out.wavelength(v) answers "what size of structure peaks at v". Reading
% out.wavelength(1) -- which looks like "the answer" and is what a caller reaches for -- returns
% the scale selected at VERTEX 1, wherever that happens to be. Planting a blob at vertex 1 then
% round-trips perfectly and the mistake is invisible: this project validated the ruler exactly
% that way (60/120/250 mm planted at vertex 1 read back 61/123/248) and carried the error into
% a tile analysis that reported a "blob size" per tile which was really the scale at one fixed,
% arbitrary vertex.
% ⭐ READ IT AT THE BLOB. Measured over 80 random seeds on this cortex, with the reading taken
% at the SEED vertex: planted 80/120/180/250 mm read 83/118/177/241, IQRs a few mm wide. Taken
% at vertex 1 the same fields read 310/310/390/431. An amplitude-weighted mean over the map is
% monotone but biased about 30% high, because the tails select coarser scales.
%
%     [~, pk] = max(abs(x));  sizeOfBlob = out.wavelength(pk);
%
% ⭐ USE THE CACHED BASIS. rheome.load.bases returns per-hemisphere Phi, Lambda and Mass already
% solved and screened -- 1000 modes down to 17 mm on this cortex -- so nothing here needs an
% eigensolver. Solving one on the spot is both slow and where the negative-eigenvalue trap
% below comes from.
%
% Filters the field with a band-pass kernel at a ladder of scales and, at every vertex, keeps
% the scale whose response is largest. That scale IS the size of the structure there: the
% standard scale-selection argument, on the surface's own operator instead of a plane.
%
% ⭐ THE NORMALISATION IS WHAT MAKES IT A SIZE. An unnormalised Laplacian response falls with
% scale, so its maximum is always at the finest scale and every blob measures as small. The
% repo's rheome.filters.mexhat is g = (t*lambda)*exp(-t*lambda), which already carries the factor t:
% it IS the scale-normalised Laplacian of Gaussian, so nothing further is needed for it. A
% custom kernel is normalised here by t^Gamma, and Gamma = 0 means "already normalised".
%
% ⚠ A VORTEX IS NOT A BLOB IN AMPLITUDE. Its detector is the winding number, which is
% invariant to any positive amplitude profile (rheome.detect.phasesingularity). Use THAT to find it
% and this to size it, by running the scale selection on the vorticity field -- not on the
% amplitude, where a strong source and a large source look the same.
%
% ⚠ AND SIZE STOPS AT THE INSTRUMENT. An inverse resolves a finite aperture; below it the
% selected scale is reading the estimator's point spread, not the cortex. MinWavelength floors
% the ladder so a size below the floor is never reported as if it were measured. On a 270-channel
% CTF array the usable floor is 122 mm, where the cortex carries about 40 independent modes.
%
% INPUTS
%   x       [nV x nT] field on the surface (real or complex; the response magnitude is used)
%   basis   .Phi [nV x K], .Lambda [K x 1], .Mass [nV x nV]   (rheome.eigen.modes)
%   Kernel  "mexhat" (default) | "diffgauss" | @(lambda, t) gains
%   Scales  the diffusion times to try; default log-spaced from Wavelengths
%   Wavelengths  [lo hi] in metres, converted to t = (w/(2*pi))^2 (default [0.03 0.5])
%   NumScales    how many, log-spaced (12)
%   MinWavelength  drop scales finer than this (0)
%   Gamma   normalisation exponent for a custom kernel (0: the kernel is already normalised)
%
% OUTPUT (struct out)
%   .scale       [nV x nT] the selected diffusion time
%   .wavelength  [nV x nT] its characteristic wavelength, 2*pi*sqrt(t), metres
%   .aperture    [nV x nT] its aperture, sqrt(2t), metres
%   .response    [nV x nT] the normalised response at the selected scale
%   .scales .wavelengths   the ladder that was searched
%   .atFloor     [nV x nT] true where the selection sat on the finest scale (unresolved)
%
% See also: rheome.detect.phasesingularity, rheome.filters.mexhat, rheome.geom.tree, rheome.eigen.modes
%
% Author: Diellor Basha, 2026

    arguments
        x
        basis (1,1) struct
        opts.Kernel = "mexhat"
        opts.Scales double = []
        opts.Wavelengths (1,2) double = [0.03 0.5]
        opts.NumScales (1,1) double {mustBeInteger, mustBePositive} = 12
        opts.MinWavelength (1,1) double {mustBeNonnegative} = 0
        opts.Gamma (1,1) double = 0
    end
    lam = double(basis.Lambda(:));
    Phi = basis.Phi;  Mass = basis.Mass;
    % ⚠ A NEGATIVE EIGENVALUE DETONATES EVERY BAND-PASS HERE. exp(-lambda*t) with lambda < 0
    % grows, so one spurious mode makes the response astronomical and the scale selection
    % pins to the coarsest rung. The cotan Laplacian can carry small negatives on obtuse
    % triangles, and eigs(...,'smallestabs') returned one of -1.6e4 on this cortex -- a solver
    % artefact, not a frequency. rheome.eigen.modes screens them; a basis from anywhere else is
    % screened here.
    bad = lam < 0;
    if any(bad)
        warning('detect:blobscale:negative', ...
            ['%d eigenvalue(s) are negative (most negative %.3g) and were clamped to zero. ' ...
             'Build the basis with rheome.eigen.modes, which screens them.'], sum(bad), min(lam));
        lam(bad) = 0;
    end
    if size(x, 1) ~= size(Phi, 1)
        error('detect:blobscale:size', 'x has %d rows; the basis is on %d vertices.', size(x,1), size(Phi,1));
    end

    t = opts.Scales(:)';
    if isempty(t)
        w = opts.Wavelengths;
        t = logspace(log10((w(1)/(2*pi))^2), log10((w(2)/(2*pi))^2), opts.NumScales);
    end
    wl = 2*pi*sqrt(t);
    keep = wl >= opts.MinWavelength;
    if ~any(keep)
        error('detect:blobscale:floor', ...
              'Every scale is finer than MinWavelength = %.0f mm; the coarsest is %.0f mm.', ...
              1e3*opts.MinWavelength, 1e3*max(wl));
    end
    t = t(keep);  wl = wl(keep);

    c = Phi' * (Mass * x);                                    % [K x nT], the modal coefficients
    nT = size(c, 2);  nV = size(Phi, 1);
    best = -inf(nV, nT);  bi = ones(nV, nT);
    for s = 1:numel(t)
        g = i_gains(opts.Kernel, lam, t(s), opts.Gamma);
        r = abs(Phi * (g .* c));                              % magnitude: complex fields too
        m = r > best;
        best(m) = r(m);  bi(m) = s;
    end
    out = struct();
    out.scale      = reshape(t(bi), nV, nT);
    out.wavelength = reshape(wl(bi), nV, nT);
    out.aperture   = sqrt(2 * out.scale);
    out.response   = best;
    out.scales     = t;
    out.wavelengths = wl;
    out.atFloor    = bi == 1;
end

% ⚠ t^Gamma, not t^(Gamma/2): the scale-normalised LAPLACIAN in two dimensions carries one
% factor of t, which is why rheome.filters.mexhat needs none added (Gamma = 0 leaves it alone).
function g = i_gains(kern, lam, t, gamma)
    if isa(kern, 'function_handle')
        g = kern(lam, t);
    else
        switch string(kern)
            case "mexhat",   g = rheome.filters.mexhat(lam, t);
            case "diffgauss", g = rheome.filters.diffgauss(lam, t, 2*t);
            otherwise
                error('detect:blobscale:kernel', ...
                      'Kernel must be "mexhat", "diffgauss" or a handle, got "%s".', string(kern));
        end
    end
    g = double(g(:));
    if gamma ~= 0, g = g * t^gamma; end
end
% Author: Diellor Basha, 2026
