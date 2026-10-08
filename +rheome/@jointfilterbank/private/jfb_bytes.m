function s = jfb_bytes(b)
% JFB_BYTES  A byte count as kB / MB / GB, for the size messages of disp and jfb_guard.
%
%   s = jfb_bytes(8.4e9)      % '8.40 GB'
%
% Author: Diellor Basha, 2026

    if     b >= 1e9, s = sprintf('%.2f GB', b/1e9);
    elseif b >= 1e6, s = sprintf('%.2f MB', b/1e6);
    else,            s = sprintf('%.1f kB', b/1e3);
    end
end

% Author: Diellor Basha, 2026
