function Y = lagfreenull(X, good, white)
% CORES.LAGFREENULL  A lag-free surrogate of the band's real signal: Gaussian noise with the band's mean spectrum, mixed by the zero-lag covariance.
%
%   Y = rheome.cores.lagfreenull(X, good, white)         % the white-noise draw is an input (reproducible)
%   Y = rheome.cores.lagfreenull(X, good, stream)        % or a RandStream: white = randn(stream, size(X))
%
% X [C x m] real band signal (the real part of rheome.cores.bandframes), good [1 x m] logical (empty: all), white [C x m]
% unit-variance Gaussian noise. amp = sqrt(mean over channels |rfft(X)|^2) per rfft bin (unscaled); R = symmetric PSD square
% root of cov(X(:,good)) (divisor n-1; eigenvalues below 0 clipped to 0); W = irfft(rfft(white) .* amp, m); each row of W is
% divided by its POPULATION std (divisor n) over the good frames; Y = single(R * W). Rows of the noise are independent, so the
% zero-lag covariance of Y matches the data's only approximately in one draw.
% nxr draws exactly one rng.standard_normal(X.shape) (numpy PCG64); a port cannot reproduce numpy's generator, so the draw is
% the input (the nxr goldens store it as null_white). With a RandStream the null is a Gaussian surrogate, not a reproducible
% quantity. Ports nxr-analysis src/nxr/vortex.py lagfree_real (line 418, commit 30b9622).
%
% See also: rheome.cores.envelopenull, rheome.cores.bandframes
%
% Author: Diellor Basha, 2026
    arguments
        X (:,:) {mustBeNumeric, mustBeReal}
        good logical
        white
    end
    X = double(X);  [C, n] = size(X);
    if isa(white, 'RandStream'), white = randn(white, C, n); end
    if ~isequal(size(white), size(X)), error('rheome:cores:lagfreenull', 'white must be %d x %d.', C, n); end
    if isempty(good), good = true(1, n); end
    good = good(:).';
    h = floor(n/2) + 1;
    amp = sqrt(mean(abs(fft(X, [], 2)).^2, 1));
    amp = amp(1:h);
    kk = 0:n-1;  full = amp(min(kk, n - kk) + 1);              % irfft: Hermitian mirror of a real gain
    A = cov(X(:, good).');                                      % rows = samples here; ddof 1
    A = (A + A.') / 2;
    [E, e] = eig(A, 'vector');
    Rh = (E .* sqrt(max(e(:).', 0))) * E.';
    W = ifft(fft(double(white), [], 2) .* full, [], 2, 'symmetric');
    W = W ./ std(W(:, good), 1, 2);
    Y = single(Rh * W);
end
