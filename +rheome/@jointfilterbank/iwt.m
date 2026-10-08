function C = iwt(obj, Y, mode)
% IWT  Reconstruct a joint spectrum from filterbank coefficients.
%
%   C = iwt(jfb, Y)             'dual' (default)
%   C = iwt(jfb, Y, 'tight')
%
% 'dual'  : wd_i = W_i / S, the canonical dual. Then sum_i conj(wd_i)*W_i = 1 wherever
%           S > 0, so iwt(wt(C)) = C EXACTLY where the bank covers. This is the one thing
%           that makes a filterbank a TRANSFORM rather than an analysis.
% 'tight' : wd_i = W_i, the frame operator. iwt(wt(C)) = S .* C.
%
% ⚠ The dual exists only where the bank COVERS. Where S is zero the dual is set to zero
% and that content is LOST -- read framebounds first. Fix the bank, not the tolerance.
%
% See also: wt, framebounds
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(mode), mode = 'dual'; end
    if ~any(strcmpi(mode, {'dual','tight'}))
        error('jointfilterbank:mode', 'mode must be ''dual'' or ''tight'', got ''%s''.', mode);
    end
    lam = obj.Lambda;  om = obj.Omega;
    K = numel(lam);  nO = numel(om);  Nf = obj.NumMembers;
    if size(Y,1) ~= K || size(Y,2) ~= nO || size(Y,3) ~= Nf
        error('jointfilterbank:size', ...
            'Y is %s but the bank is [%d %d %d].', mat2str(size(Y)), K, nO, Nf);
    end

    useDual = strcmpi(mode, 'dual');
    if useDual
        S    = jfb_frameop(obj);
        bad  = S <= 1e-12 * max([max(S(:)); realmin]);
        Sinv = zeros(size(S));
        Sinv(~bad) = 1 ./ S(~bad);
    end

    C = zeros(K, nO);
    if ~isreal(Y), C = complex(C); end
    for m = 1:Nf
        Wm = jfb_member(obj, m, lam, om);
        if useDual, Wm = Wm .* Sinv; end
        C = C + conj(Wm) .* Y(:,:,m);
    end
end

% Author: Diellor Basha, 2026
