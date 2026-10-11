function h = hash(varargin)
% MODES.HASH  SHA-256 (hex) of numeric arrays: the fingerprint of frozen centroids.
%
%   h = rheome.modes.hash(C1, C2, ...)       % complex arrays hashed as [real; imag], double, in order
%
% Size and class do not change the hash of equal values: each array enters as double [real(:); imag(:)],
% preceded by its size, so a reshape is a different hash.
%
% Author: Diellor Basha, 2026

    md = java.security.MessageDigest.getInstance('SHA-256');
    for k = 1:numel(varargin)
        x = varargin{k};
        md.update(typecast(double(size(x)), 'uint8'));
        md.update(typecast([real(double(x(:))); imag(double(x(:)))], 'uint8'));
    end
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, []));
end

% Author: Diellor Basha, 2026
