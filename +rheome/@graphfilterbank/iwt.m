function X = iwt(obj, W, mode)
% IWT  Inverse graph wavelet transform.
%
%   X = iwt(gfb, W)            'dual'  (default)
%   X = iwt(gfb, W, 'tight')
%
% 'dual'  : w_m = g_m / S, the canonical dual. iwt(wt(X)) = X EXACTLY for ANY family
%           wherever the bank covers -- unlike icwt, which needs an analytic wavelet.
%           This is what makes a non-tight scale-space an invertible transform.
% 'tight' : w_m = g_m, the frame operator. iwt(wt(X)) = A*X, exact only when tight.
%
% See also: wt, framebounds, isframetight
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(mode), mode = 'dual'; end
    if ~any(strcmpi(mode, {'dual','tight'}))
        error('graphfilterbank:mode', 'mode must be ''dual'' or ''tight'', got ''%s''.', mode);
    end
    T   = gfb_transform(obj);
    lam = gfb_lambda(obj, T);
    H   = graphfilters(obj, 'Lambda', lam, 'Dual', strcmpi(mode,'dual'));

    M = obj.NumMembers;  nT = size(W, 2);
    X = zeros(size(W,1), nT);
    if ~isreal(W), X = complex(X); end
    for m = 1:M
        X = X + T.inverse(H(:, m) .* T.forward(W(:, :, m)));
    end
end

% Author: Diellor Basha, 2026
