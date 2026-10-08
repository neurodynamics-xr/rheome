function b = framebounds(obj, lambda)
% FRAMEBOUNDS  Frame bounds A, B and tightness. No cwtfilterbank analogue.
%
%   b = framebounds(gfb)           over the spectral interval (a dense grid)
%   b = framebounds(gfb, Lambda)   EXACT, on the spectrum the transform will see
%
% A = B is TIGHT: synthesis(analysis(F)) = A*F. A -> 0 means part of the spectrum is
% uncovered and the canonical dual does not exist there.
%
% Bounds over the continuum are NOT the bounds the transform sees -- pass the spectrum
% when you have one.
%
% OUTPUT: b.A  b.B  b.Tightness (B/A)  b.Uncovered (fraction with S ~ 0)
%
% See also: isframetight, graphfilters
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(lambda)
        if obj.HasSpectrum, lambda = obj.Lambda_;
        else,               lambda = linspace(0, obj.Lmax_, 512)';
        end
    end
    H = graphfilters(obj, 'Lambda', lambda);
    S = sum(H.^2, 2);
    b = struct('A', min(S), 'B', max(S), ...
               'Tightness', max(S) / max(min(S), eps), ...
               'Uncovered', mean(S <= eps * max([max(S); realmin])));
end

% Author: Diellor Basha, 2026
