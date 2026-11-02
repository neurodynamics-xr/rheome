function y = boxsmooth(x, k)
% CORES.BOXSMOOTH  Moving mean with zero-padded edges, numpy's np.convolve(x, ones(k)/k, 'same').
%
%   y = rheome.cores.boxsmooth(x, k)
%
% Not movmean: movmean shrinks its window at the ends, this divides by k there and counts the missing samples as 0.
% For even k the window of output i is x[i-k/2 ... i+k/2-1] (one sample earlier than centred); odd k is symmetric.
% MATLAB's conv(...,'same') would start one sample later for even k, so the central part is cut here as numpy cuts it.
% x is a vector (any orientation, kept); k a positive integer no longer than x.
%
% Ports the smoothing of nxr-analysis src/nxr/vortex.py burst_stages, spindle_envelope, risefall and envelope_matched_null
% (commit 30b9622).
%
% See also: rheome.cores.burstenvelope, rheome.cores.burstfit
%
% Author: Diellor Basha, 2026
    arguments
        x {mustBeNumeric, mustBeVector}
        k (1,1) double {mustBeInteger, mustBePositive}
    end
    n = numel(x);
    if k > n, error('rheome:cores:boxsmooth', 'window k = %d is longer than the signal (%d).', k, n); end
    full = conv(double(x(:)).', ones(1, k) / k);
    y = full((1:n) + floor((k - 1) / 2));
    if iscolumn(x), y = y(:); end
end
