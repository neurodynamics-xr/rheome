function G = travwave(lambda, f, c, width, sense, half)
% FILTERS.TRAVWAVE  Joint traveling-wave filter: energy on the dispersion ridge, with a sense.
%
%   G = rheome.filters.travwave(lambda, f, c, width)              % both senses (|f|)
%   G = rheome.filters.travwave(lambda, f, c, width, sense)       % sense = +1 | -1 | 0
%   G = rheome.filters.travwave(lambda, f, c, width, sense, half) % half = 'full' (default) | 'positive'
%
% A non-separable joint filter g(lambda,f) selecting the (lambda,f) pairs that obey a linear
% dispersion relation -- waves propagating on the cortex at phase speed c. Energy concentrates
% on the ridge  f = sense * c*sqrt(lambda)/(2*pi)  with a Gaussian width.
%
% UNITS. f is the SIGNED temporal frequency in HERTZ; lambda in 1/length^2 and c in length/s, so
% c*sqrt(lambda)/(2*pi) is in hertz; width in hertz.
%
% SENSE. The sign of f is the sense of travel of a complex (rotary) coefficient, so a wave travelling
% one way and the same wave travelling the other way sit on opposite sides of f = 0:
%   sense = +1  the ridge at +c*sqrt(lambda)/(2*pi) only (positive frequencies)
%   sense = -1  the ridge at -c*sqrt(lambda)/(2*pi) only
%   sense =  0  both ridges, i.e. a function of |f| (default; the right choice for the real
%               coefficients of a real field, whose spectrum is even in f)
% A sense needs a two-sided frequency axis: on a one-sided ('positive') grid only sense = +1 and 0
% have bins to act on.
%
% INPUTS:
%   lambda [K x 1] eigenvalues (dbasis.Lambda)
%   f      [1 x N] frequency axis (Hz); a full FFT grid (rheome.filters.jspectrum) unless half says otherwise
%   c      phase speed (length/s), default 1
%   width  ridge width (Hz), default 2
%   sense  +1 | -1 | 0, default 0
%   half   'full' (default) | 'positive': see rheome.filters.signedfrequencies
%
% OUTPUT:
%   G      [K x N] joint gains (multiply the joint spectrum elementwise)
%
% See also: rheome.filters.signedfrequencies, rheome.filters.jspectrum, rheome.filters.stmatern, rheome.filters.bandpass
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(c),     c = 1;     end
    if nargin < 4 || isempty(width), width = 2; end
    if nargin < 5 || isempty(sense), sense = 0; end
    if nargin < 6 || isempty(half),  half = 'full'; end
    if ~any(sense == [-1 0 1])
        error('filters:travwave:sense', 'sense must be +1, -1 or 0, got %g.', sense);
    end
    ws = rheome.filters.signedfrequencies(f, half);
    fr = c .* sqrt(max(double(lambda(:)), 0)) ./ (2*pi);   % [K x 1] ridge frequency, Hz
    if sense == 0
        fx = abs(ws);   G = exp( -((fx - fr).^2) / (2 * max(width,eps)^2) );
    else
        G = exp( -((ws - sense*fr).^2) / (2 * max(width,eps)^2) );
    end
end

% Author: Diellor Basha, 2026
