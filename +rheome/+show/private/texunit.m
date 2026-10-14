function s = texunit(u)
% TEXUNIT  A unit as LaTeX maths (tokens type.math: upright, exponents in maths): 'nA·m' -> '$\mathrm{nA{\cdot}m}$',
% 'fT^2' -> '$\mathrm{fT^{2}}$', 'µV' -> '$\mathrm{\mu V}$'. A unit already holding $ is returned as given.
%
% Author: Diellor Basha, 2026
    s = char(u);
    if contains(s, '$') || isempty(s) || all(s < 128 & s ~= '^'), return; end   % plain ASCII units stay text
    s = regexprep(s, '\^(-?\d+)', '^{$1}');
    s = strrep(s, char(178), '^{2}');  s = strrep(s, char(179), '^{3}');
    s = strrep(s, char(183), '{\cdot}');  s = strrep(s, char([181]), '\mu ');  s = strrep(s, char(956), '\mu ');
    s = strrep(s, '%', '\%');  s = strrep(s, ' ', '\,');
    s = ['$\mathrm{' s '}$'];
end
