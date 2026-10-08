function W = jointfilters(obj, m, varargin)
% JOINTFILTERS  The bank's filters on the joint plane.
%
%   W = jointfilters(jfb, m)              one member  [K x nOmega]      -- always safe
%   W = jointfilters(jfb)                 ALL members [K x nOmega x Nf] -- GUARDED
%   W = jointfilters(jfb, [], 'Force', true)
%
% ⭐ The bank stores factors, not products. A single member is cheap and always
% available; materialising all of them is an explicit request and refuses above
% MaxBytes, because [K x nOmega x Nf] is 8.4 GB for a broadband cortex analysis with
% speed selectivity.
%
% See also: framebounds, wt, scalogram, spectralResponse
%
% Author: Diellor Basha, 2026

    if nargin >= 2 && ~isempty(m)
        W = jfb_member(obj, m, obj.Lambda, obj.Omega);
        return;
    end

    p = inputParser;
    p.addParameter('Force', false);
    p.parse(varargin{:});

    lam = obj.Lambda;  om = obj.Omega;
    K = numel(lam);  nO = numel(om);  Nf = obj.NumMembers;
    if ~p.Results.Force
        jfb_guard(obj, [K nO Nf], 'the whole bank');
    end
    W = zeros(K, nO, Nf);
    for mm = 1:Nf
        W(:,:,mm) = jfb_member(obj, mm, lam, om);
    end
end

% Author: Diellor Basha, 2026
