function d = index(n, opts)
% DOMAIN.INDEX  A one-dimensional index domain: time, frequency, or an eigenmode axis.
%
%   d = rheome.domain.index(360000, Name="time", Units="s")
%   d = rheome.domain.index(400, Name="lb_modes", Of=dCortex)
%
% ⭐ THE SPECTRAL LINE IS A DOMAIN TOO. Eigenmode coefficients are a field on the mode
% index, not on the cortex -- which is the distinction that lets the flow pipeline work
% entirely in coefficients and lift to vertices only to draw. `Of` records which domain the
% line was derived from, so a coefficient field still knows the surface behind it.
%
% See also: rheome.domain.of, rheome.domain.kinds
%
% Author: Diellor Basha, 2026

    arguments
        n (1,1) double {mustBeInteger, mustBePositive}
        opts.Name  (1,1) string = ""
        opts.Units (1,1) string = ""
        opts.Of    = []
    end
    d = struct('kind', 'index', 'nV', n, 'nE', 0, 'nF', 0, 'id', '', ...
               'name', char(opts.Name), 'source', '', 'units', char(opts.Units), 'of', '');
    if ~isempty(opts.Of)
        o = rheome.domain.of(opts.Of);
        d.of = o.id;
    end
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(sprintf('index|%d|%s|%s', n, d.name, d.of)));
    dg = typecast(md.digest(), 'uint8');
    d.id = lower(reshape(dec2hex(dg(1:4), 2)', 1, []));
end
% Author: Diellor Basha, 2026
