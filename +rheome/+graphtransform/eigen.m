function T = eigen(Phi, B, Lambda, varargin)
% GRAPHTRANSFORM.EIGEN  Exact-eigenbasis transform for a graphfilterbank.
%
%   T = rheome.graphtransform.eigen(Phi, B, Lambda)
%   T = rheome.graphtransform.eigen(Phi, B, Lambda, 'nV', nV)   % quaternion: 4*nV rows
%
% Everything operator-specific lives behind four handles:
%   .forward(F) = Phi' * (B * F)      .inverse(C) = Phi * C
%   .norm(F)    per-vertex magnitude  .lambda     the exact spectrum
%
% For a quaternion (Dirac) basis pass 'nV' and .norm measures the IMAGINARY (x,y,z)
% rows only -- physical current magnitude, not the full quaternion norm. That is the
% only thing Dirac changes; the filtering itself is identical.
%
% ⭐ THE DIRAC SPECTRUM IS 4-FOLD DEGENERATE, AND THE FOUR ARE ONE PATTERN. Measured on the
% cached a reference subject basis: eigenvalues arrive in groups of four agreeing to 4.8e-14 relative,
% with consecutive groups 2.1% apart. Those four are a QUATERNIONIC LINE -- psi*i, psi*j and
% psi*k all lie in the span of the group to a residual of 1e-13 -- because the operator is
% right-H-linear, so every eigenspace is an H-line, four REAL dimensions carrying one spatial
% scale. The solver returns an arbitrary real-orthonormal basis of it, so the four vectors you
% get are four arbitrary quaternion rotations of the same pattern and are not separately
% meaningful. ⚠ Count SCALES in groups of four, never modes.
%
% ⚠⚠ THE FOUR COEFFICIENTS WITHIN A GROUP ARE A GAUGE, NOT FOUR MEASUREMENTS. It is tempting to
% read them as the w/x/y/z parts of the field at that scale and to recover direction from them.
% They are not, and you cannot. Measured: the cached basis is NOT stored in canonical
% [psi, psi*i, psi*j, psi*k] order -- ||Q2 - psi*i|| comes out at 1.4-1.7, and sqrt(2) is what
% two random orthonormal vectors give, so the stored columns are an arbitrary orthonormal basis
% of the H-line. Remixing one group by an arbitrary orthogonal R changes the individual
% coefficients by 122% of the group norm while leaving |c| invariant exactly and the
% RECONSTRUCTION invariant to 2.8e-16. Any orthonormal basis of a degenerate eigenspace is a
% legitimate solver output, so within a group the ONLY invariant is the total power.
% ⭐ DIRECTION LIVES IN THE RECONSTRUCTION, and it is exact there: filtering a seeded current of
% [1 0.4 0] and reading the field at the seed returns [1.004 0.391 0.004] under one gauge and
% [1.004 0.390 0.004] under another. Reconstruct at each scale and read the 3-vector per vertex;
% that gives direction AS A FUNCTION OF SCALE, which is what the coefficients cannot.
% ⚠ If per-component meaning is wanted, the gauge has to be IMPOSED -- rotate each group onto a
% chosen reference -- and it is then a property of that choice, not of the basis.
%
% ⚠ AND A FILTERED PURE-IMAGINARY FIELD DOES NOT STAY PURE-IMAGINARY, which is why .norm
% discards something. Seeding a delta with a pure current 3-vector and applying a mexhat in this
% basis leaks energy into the REAL (w) row: measured 1.2% at sigma 250 mm, 4.5% at 150 mm and
% 9.7% at 100 mm -- it grows as the filter sharpens. So a "Dirac wavelet" is a 3-vector field to
% 90-99%, not exactly. ⭐ The leak is not at the seed: seeding [1 0 0] and reading the wavelet AT
% the seed vertex returns [1.000 -0.011 0.003] with w = -0.000, so the direction is preserved
% where it was planted and the scalar part lives in the lobes.
%
% Phi must be B-orthonormal (Phi'*B*Phi = I); Parseval depends on it, exactly as
% cwtfilterbank depends on fft's.
%
% See also: rheome.graphtransform.chebyshev, rheome.graphtransform.validate, rheome.graphfilterbank
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('nV', []);
    p.parse(varargin{:});
    nV = p.Results.nV;

    rows = size(Phi, 1);
    isQuat = ~isempty(nV) && rows == 4*nV;
    if isempty(nV), nV = rows; end

    if isQuat
        iImag = reshape((0:nV-1)*4 + (2:4)', [], 1);     % (x,y,z) rows = physical current
        nrm = @(F) i_quatnorm(F, iImag, nV);
    else
        nrm = @(F) abs(F);
    end

    T = struct('forward', @(F) Phi' * (B * F), ...
               'inverse', @(C) Phi * C, ...
               'norm',    nrm, ...
               'lambda',  double(Lambda(:)), ...
               'nV',      nV, ...
               'rows',    rows);
end

function e = i_quatnorm(F, iImag, nV)
    G = reshape(F(iImag, :), 3, nV, []);
    e = reshape(sqrt(sum(abs(G).^2, 1)), nV, []);
end

% Author: Diellor Basha, 2026
