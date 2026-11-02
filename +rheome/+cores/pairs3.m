function s = pairs3(fa, fb, Sx)
% CORES.PAIRS3  Poincare S3 of the pair (fa, fb) under Sx: 2 Im(fa Sx fb') / Re(fa Sx fa' + fb Sx fb'), in -1 ... 1 (nxr s3).
%
%   s = rheome.cores.pairs3(fa, fb, Sx)
%
% fa, fb real row (or column) vectors of length C; invariant to rotations in the (fa, fb) plane. Uses the plain transposes of the
% real filters, as nxr (fa @ Sx @ fb).
%
% Ports nxr-analysis src/nxr/vortex.py s3 (:148, commit 3da1320).
%
% See also: rheome.cores.pair, rheome.cores.pairshare
%
% Author: Diellor Basha, 2026
    arguments
        fa (:,:) double
        fb (:,:) double
        Sx (:,:) {mustBeNumeric}
    end
    fa = fa(:).'; fb = fb(:).';
    s = 2 * imag(fa * Sx * fb.') / real(fa * Sx * fa.' + fb * Sx * fb.');
end
