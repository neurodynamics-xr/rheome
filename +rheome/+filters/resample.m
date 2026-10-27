function [y, p, q] = resample(x, fs, to, opts)
% FILTERS.RESAMPLE  Rows of a recording resampled from fs to another rate, anti-aliased (polyphase, a Kaiser-windowed sinc).
%
%   y = rheome.filters.resample(x, fs, to)            % x [C x T] at fs Hz -> y [C x ceil(T*p/q)] at `to` Hz
%   [y, p, q] = rheome.filters.resample(x, 1200, 300, Dim=2)
%
% The ratio to/fs is the fraction p/q of the two rates in mHz, reduced (both at most 1000: a rate pair with no small ratio is
% refused). The samples are upsampled by p, filtered by a low-pass FIR and downsampled by q in one polyphase pass (upfirdn);
% the filter is the ideal low-pass at 1/max(p, q) of Nyquist, 2*10*max(p, q) + 1 taps, Kaiser window beta = 5, scaled to unit
% gain at DC and by p; the output is trimmed by the filter's delay so sample j of y is at time j/to, as sample i of x is at i/fs.
% Beyond the input's ends the signal is taken as zero (scipy's padtype 'constant'): the first and last ~10*max(p,q)/q output
% samples see the edge, so a chunk is resampled with context and trimmed, or the whole recording at once.
%
% Port of scipy.signal.resample_poly (window ('kaiser', 5.0), padtype 'constant'), what nxr-convert resample_rows uses: the
% 300 Hz copy of a stored recording is this function's output. MATLAB's own resample designs its filter differently (firls).
%
% Example:
%   y = rheome.filters.resample(R.Data, R.Fs, 300);     % R = rheome.io.read.nxrrecording(...)
%
% See also: rheome.io.read.nxrrecording, rheome.spectral.root, upfirdn
%
% Author: Diellor Basha, 2026 (port of scipy.signal.resample_poly)
arguments
    x {mustBeNumeric}
    fs (1,1) double {mustBePositive}
    to (1,1) double {mustBePositive}
    opts.Dim (1,1) double {mustBeMember(opts.Dim, [1 2])} = 2
end
p = round(to * 1000); q = round(fs * 1000); g = gcd(p, q); p = p / g; q = q / g;
if p > 1000 || q > 1000, error('rheome:resample', 'cannot resample %g Hz to %g Hz: the ratio %d/%d is not a small fraction', fs, to, p, q); end
if opts.Dim == 2, x = x.'; end                                  % upfirdn works down columns
cls = class(x); x = double(x);
if p == q, y = x; else, y = i_poly(x, p, q); end
y = cast(y, cls);
if opts.Dim == 2, y = y.'; end
end

function y = i_poly(x, p, q)
nIn = size(x, 1); nOut = ceil(nIn * p / q);
m = max(p, q); half = 10 * m; L = 2 * half + 1;
n = (0:L-1)' - half;
h = (1 / m) * i_sinc(n / m) .* kaiser(L, 5);                    % firwin(L, 1/m, window=('kaiser', 5.0))
h = h / sum(h) * p;
pre = q - mod(half, q); drop = floor((half + pre) / q);
post = 0;
while i_outlen(L + pre + post, nIn, p, q) < nOut + drop, post = post + 1; end
h = [zeros(pre, 1); h; zeros(post, 1)];
y = upfirdn(x, h, p, q);
y = y(drop + (1:nOut), :);
end

function n = i_outlen(nh, nIn, p, q)
% MATLAB upfirdn's output length (trailing zeros in h only lengthen it: they change no value kept)
n = ceil(((nIn - 1) * p + nh) / q);
end

function s = i_sinc(x)
s = ones(size(x)); k = x ~= 0; s(k) = sin(pi * x(k)) ./ (pi * x(k));
end
