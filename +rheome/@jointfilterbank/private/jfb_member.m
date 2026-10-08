function W = jfb_member(obj, m, lambda, omega)
% JFB_MEMBER  Evaluate ONE member on the joint grid -> [K x nOmega].
%
%   W = jfb_member(obj, m, lambda, omega)
%
% ⭐ The kernel is called as K(lambda [K x 1], omega [1 x nOmega]) and relies on IMPLICIT
% EXPANSION. Never build ndgrid(lambda, omega): two full [K x nOmega] grids is 338 MB at
% cortex scale, which defeats the whole point of storing factors.
%
% Author: Diellor Basha, 2026

    idx = obj.Index;
    if ~isscalar(m) || m < 1 || m > size(idx,1) || mod(m,1) ~= 0
        error('jointfilterbank:member', ...
            'member must be an integer in 1..%d, got %s.', size(idx,1), mat2str(m));
    end
    ig = idx(m,1);  it = idx(m,2);  ik = idx(m,3);

    lambda = double(lambda(:));
    omega  = double(omega(:)).';
    K = numel(lambda);  nO = numel(omega);

    gv = obj.GraphBank.gain(ig);  gv = gv(lambda);  gv = gv(:);        % [K x 1]
    hv = obj.TimeMembers{it}(omega);  hv = hv(:).';                    % [1 x nOmega]
    Kv = obj.KernelMembers{ik}(lambda, omega);

    if isscalar(Kv)
        W = Kv * (gv * hv);
    else
        if ~isequal(size(Kv), [K nO])
            error('jointfilterbank:kernelSize', ...
                ['Joint kernel %d returned %s but the grid is %s. A kernel is called as ' ...
                 'K(lambda [%d x 1], omega [1 x %d]) and must return a scalar or that ' ...
                 'size by implicit expansion.'], ik, mat2str(size(Kv)), ...
                mat2str([K nO]), K, nO);
        end
        W = Kv .* (gv * hv);
    end
end

% Author: Diellor Basha, 2026
