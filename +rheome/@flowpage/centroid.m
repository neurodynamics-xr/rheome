function k = centroid(obj, iFrame, basis)
% CENTROID  Energy-weighted mean wavenumber at one frame -- the spatial scale, in rad/m.
%
%   k = centroid(fp, iFrame)         2*pi/k is the scale in metres
%
% ⭐ AMPLITUDE-INVARIANT BY CONSTRUCTION. Scaling the field scales every scale's energy
% equally and leaves the centroid unchanged, so any dependence on band power is STRUCTURAL
% rather than mechanical. That is what makes it usable as a descriptor against power.
%
% ⚠ TWO BASES, AND THE DEFAULT MATTERS. 'mode' (default) weights sqrt(lambda) by the raw
% mode energy; 'bank' weights the graph members' centre wavenumbers by their energies.
%
% MEASURED on 4000 core frames of resting alpha: the BANK centroid moved over only 55-58 mm
% while the mode centroid moves freely -- with 7 heavily overlapping members the bank
% centroid is pinned near the bank's own centre of mass and reports the FILTERBANK, not the
% data. Use 'bank' only where a hard ROI mask would otherwise be needed: masking a field by
% a region indicator is a rectangular spatial window whose spectral leakage reads as
% fine-scale structure, which the bank's smooth windows avoid. For a WHOLE-CORTEX centroid
% there is no mask, so 'mode' is both correct and better resolved.
%
% See also: scaleEnergy, globalIndex
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(basis), basis = 'mode'; end
    switch lower(basis)
        case 'mode'
            E  = modeSpectrum(obj, iFrame).';
            kc = sqrt(obj.Lambda(:)).';
        case 'bank'
            E  = scaleEnergy(obj, iFrame);
            kc = double(reshape(centerWavenumbers(obj.GraphBank), 1, []));
        otherwise
            error('flowpage:basis', 'basis must be ''mode'' or ''bank'', got ''%s''.', basis);
    end
    w = sum(E);
    if w <= 0, k = NaN; return; end
    k = sum(kc .* E) / w;
end

% Author: Diellor Basha, 2026
