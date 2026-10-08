function X = read(obj, a, b)
% READ  Raw partial read of samples a..b for the selected channels.
%
%   X = read(pr, a, b)      -> [NumChannels x (b-a+1)] in pr.Precision
%
% The columns come off disk as a hyperslab; the channel subset is taken in memory
% afterwards, because F is column-major and a few hundred rows cost nothing next to a
% strided HDF5 read.
%
% See also: page, pagerange
%
% Author: Diellor Basha, 2026

    if a < 1 || b > obj.NumSamples || b < a
        error('pagedrecording:range', ...
            'Sample range [%d %d] is outside 1..%d.', a, b, obj.NumSamples);
    end

    X = obj.MF_.F(:, a:b);
    X = X(obj.Channels, :);

    if strcmp(obj.Precision, 'single'), X = single(X); else, X = double(X); end
end

% Author: Diellor Basha, 2026
