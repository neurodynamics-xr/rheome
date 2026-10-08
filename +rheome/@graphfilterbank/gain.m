function h = gain(obj, m)
% GAIN  The gain handle g_m(lambda) of one member.
%   h = gain(gfb, m)      m in 1..NumMembers
% Author: Diellor Basha, 2026
    if ~isscalar(m) || m < 1 || m > obj.NumMembers || mod(m,1) ~= 0
        error('graphfilterbank:member', ...
            'member must be an integer in 1..%d, got %s.', obj.NumMembers, mat2str(m));
    end
    h = obj.G_{m};
end

% Author: Diellor Basha, 2026
