function n = tfb_smooth(n)
% TFB_SMOOTH  Smallest integer >= n whose prime factors are all <= 5.
%
% Sub-band FFT lengths are rounded up to a 5-smooth number because MATLAB's FFT is fastest
% there, and that rounding is what fixes each member's rate. `wt` and `rates` must round
% identically, or a member's instants would disagree with its declared rate -- hence one
% definition, here.
%
% Author: Diellor Basha, 2026
    while true
        r = n;
        for p = [2 3 5], while mod(r, p) == 0, r = r / p; end, end
        if r == 1, return; end
        n = n + 1;
    end
end
% Author: Diellor Basha, 2026
