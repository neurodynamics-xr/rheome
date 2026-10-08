function coeff = mode_coefficients(datasetName, basis, opts)
% DYNAMICS.MODE_COEFFICIENTS  Reduced-order coefficient time series for the dynamics fits.
%
%   coeff = rheome.dynamics.mode_coefficients(datasetName, basis [,opts])
%
% Projects the reconstructed cortical activity onto a spatial eigenbasis, giving the compact
% coefficient series C [K x T] that rheome.dynamics.dmd / dispersion / pde_fit operate on.
%   basis = 'lbo'   : Laplace–Beltrami coefficients of the source magnitude |J| (clean
%                     spatial-scale axis; the interpretable, dispersion-friendly choice).
%   basis = 'dirac' : the Dirac vector mode-coefficient series out.c (full vector dynamics).
%
% OUTPUT (struct coeff):
%   .C [K x T]           coefficient series
%   .dt                  sampling interval (s)
%   .lambda [K x 1]      per-mode spatial eigenvalue (LBO) -- the dispersion x-axis; for 'dirac' this
%                        is the scalar-LBO reference lambdaLBO (Dirac lambda is co-normalised)
%   .basisModes          [nV x K] spatial modes (for reconstruction to the cortex)
%
% See also: rheome.source.dirac, rheome.eigen.modes, rheome.dynamics.dmd
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'K') || isempty(opts.K), opts.K = 400; end

    dataset = rheome.load.dataset(datasetName);
    coeff.dt = 1 / dataset.rec.sfreq;

    switch lower(basis)
        case 'dirac'
            out = rheome.source.dirac(datasetName, [], [], opts.K, 'amplitude');
            coeff.C          = out.c;                       % [nModes x T]
            coeff.basisModes = out.dbasis.Phi;
            % per-Dirac-mode eigenvalue (dimension matches C). NOTE co-normalised (scale-invariant,
            % [0,1]) -- the dispersion axis is therefore co-normalised for the Dirac basis; the LBO
            % basis gives the physical (1/m) scale.
            coeff.lambda     = out.dbasis.Lambda(:);
        case 'lbo'
            current = dataset.inverse.ImagingKernel * dataset.rec.F(dataset.inverse.GoodChannel, :);
            magnitude = sqrt(current(1:3:end,:).^2 + current(2:3:end,:).^2 + current(3:3:end,:).^2);
            coeff.C          = dataset.basis.Phi' * (dataset.basis.Mass * magnitude);
            coeff.basisModes = dataset.basis.Phi;
            coeff.lambda     = dataset.basis.Lambda(:);
        otherwise
            error('dynamics:mode_coefficients:basis', 'basis must be ''lbo'' or ''dirac''.');
    end
end

% Author: Diellor Basha, 2026
