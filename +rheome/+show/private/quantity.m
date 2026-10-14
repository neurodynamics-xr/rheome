function q = quantity(q, x)
% QUANTITY  The quantity whose colormap a field takes (tokens quantity: one colormap per quantity). Given
% q is checked against the tokens; empty q is inferred from the data x: signed -> 'divergence' (the
% tokens' alias of the retired 'signed'), one-signed -> 'magnitude'. Name it to get its own map.
%
% Author: Diellor Basha, 2026
    T = tokens();
    if isempty(q)
        if any(x(:) < 0) && any(x(:) > 0), q = 'divergence'; else, q = 'magnitude'; end
    end
    q = char(q);
    if ~isfield(T.quantity, q) || strcmp(q, 'aliases') || strcmp(q, 'never') || strcmp(q, 'anatomy')
        error('rheome:show:quantity', 'no colormap for quantity ''%s'' (tokens quantity); ask Graphics for one', q);
    end
end
