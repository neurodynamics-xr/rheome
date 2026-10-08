function hAx = spectrum2d(lambda, freqs, P, varargin)
% SHOW.SPECTRUM2D  Joint eigenmode-frequency map: power (or gain) over (sqrt(lambda), Hz).
%
%   rheome.show.spectrum2d(lambda, freqs, P)
%   rheome.show.spectrum2d(lambda, freqs, P, name, value, ...)
%
% Displays a joint time-vertex quantity as a 2-D image: temporal frequency (Hz) on x,
% spatial frequency sqrt(lambda) on y (row = mode, ascending lambda). Use for the joint
% spectrum |chat|^2 or a joint-filter gain g(lambda,omega).
%
% INPUTS:
%   lambda [K x 1] eigenvalues (one per row of P; ascending)
%   freqs  [1 x N] one-sided frequency axis (Hz) from rheome.filters.jspectrum
%   P      [K x N] the quantity to show (power, or gain)
%
% OPTIONS:
%   'FreqRange' [lo hi] Hz window on x (default [0 40])
%   'Log'       true -> log10 color scale (default true; good for power, off for gains)
%   'Colormap'  default parula
%   'Parent','Title','Visible'
%
% OUTPUT: hAx axes handle.
%
% See also: rheome.filters.jspectrum, rheome.filters.travwave, rheome.filters.stmatern
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('FreqRange', [0 40]);
    p.addParameter('Log', true);
    p.addParameter('Colormap', parula(256));
    p.addParameter('Parent', []);
    p.addParameter('Title', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    opt = p.Results;

    lambda = lambda(:);
    keep = (freqs >= opt.FreqRange(1)) & (freqs <= opt.FreqRange(2));
    Pk = P(:, keep);  fk = freqs(keep);
    if opt.Log, Pk = log10(Pk + eps); end

    if isempty(opt.Parent)
        hFig = figure('Color','w','Visible',opt.Visible);
        hAx  = axes('Parent', hFig);
    else
        hAx = opt.Parent;
    end

    K = numel(lambda);
    imagesc(hAx, fk, 1:K, Pk);
    set(hAx, 'YDir', 'normal');
    colormap(hAx, opt.Colormap);  colorbar(hAx);
    xlabel(hAx, 'temporal frequency (Hz)');
    ylabel(hAx, 'spatial frequency  \surd\lambda  (1/m)');

    % relabel the y ticks (mode index -> sqrt(lambda), the physical spatial frequency)
    yt = unique(round(linspace(1, K, 6)));
    set(hAx, 'YTick', yt, 'YTickLabel', arrayfun(@(k) sprintf('%.0f', sqrt(max(lambda(k),0))), yt, 'uni', 0));
    if ~isempty(opt.Title), title(hAx, opt.Title, 'Interpreter','none'); end

    if nargout == 0, clear hAx; end
end

% Author: Diellor Basha, 2026
