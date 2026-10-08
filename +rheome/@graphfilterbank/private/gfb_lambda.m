function lam = gfb_lambda(obj, T)
% GFB_LAMBDA  The spectrum the gains should be evaluated on.
% Author: Diellor Basha, 2026
    if isfield(T,'lambda') && ~isempty(T.lambda)
        lam = T.lambda;
    elseif obj.HasSpectrum
        lam = obj.Lambda_;
    else
        error('graphfilterbank:noSpectrum', ...
            ['The Transform carries no .lambda and the bank was not built from a ' ...
             'spectrum, so there is nothing to evaluate the gains on.']);
    end
end

% Author: Diellor Basha, 2026
