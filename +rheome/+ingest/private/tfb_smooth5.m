function n = tfb_smooth5(n)
% TFB_SMOOTH5  Smallest integer >= n whose prime factors are all <= 5 (the FFT length rule).
% The same rounding @timefilterbank uses, repeated here because a package cannot reach into
% another class's private folder; the two must agree or a pair's grid would not be a member's.
% Author: Diellor Basha, 2026
    while true
        r = n;
        for p = [2 3 5], while mod(r, p) == 0, r = r / p; end, end
        if r == 1, return; end
        n = n + 1;
    end
end
% Author: Diellor Basha, 2026
