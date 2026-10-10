function d = designbank(family, Nf, tmin, tmax, lmax, lminEff, voices)
% DESIGNBANK  Build the member gain handles for one family.
%   d = designbank(family, Nf, tmin, tmax, lmax, lminEff, voices)
%
% Lifted verbatim from rheome.filters.frame so behaviour is preserved bit-for-bit. The scale
% RANGE arrives already resolved (tmin/tmax); Nf controls only DENSITY.
%
% OUTPUT: d.g {1 x M} gain handles, d.t [1 x M] scale param (NaN where none), d.Nf
%
% Author: Diellor Basha, 2026

    switch lower(family)
        case 'itersine'
            if isempty(Nf), Nf = 8; end
            overlap = 2;
            scale   = lmax / (Nf - overlap + 1) * overlap;
            kf = @(x) sin(0.5*pi*(cos(pi*x)).^2) .* (x >= -0.5 & x <= 0.5);
            g = cell(1, Nf);
            for ii = 1:Nf
                g{ii} = @(l) kf(double(l(:))/scale - (ii-overlap/2)/overlap) ./ sqrt(overlap) .* sqrt(2);
            end
            t = nan(1, Nf);

        case 'logitersine'
            % Itersine windows on the LOG-WAVENUMBER axis xi = voices*log2(sqrt(lambda)/k_top),
            % two voices wide, member j centred at xi = -j, plus a low-pass and a high-pass
            % member completing the frame: sum_m g_m^2 = 1 on [0, lmax]. Tight AND
            % log-spaced -- the temporal timefilterbank's construction on a graph spectrum.
            % Bands in octaves of wavelength.
            V = voices;  kTop = sqrt(lmax);  kMin = sqrt(max(lminEff, eps));
            J = max(1, floor(V * log2(kTop / kMin)));            % lowest centre >= xi(kMin)
            xi = @(l) V * log2(max(sqrt(double(l(:))), realmin) / kTop);
            w  = @(u) sin(0.5*pi * cos(0.5*pi * u).^2) .* (abs(u) <= 1);
            g = cell(1, J + 2);
            g{1} = @(l) sqrt(max(1 - w(xi(l) + J).^2, 0)) .* (xi(l) < -J);   % 1 at lambda = 0 (xi -> -inf)
            for jj = 1:J
                g{J + 2 - jj} = @(l) w(xi(l) + jj);              % ascending in k: coarse first
            end
            g{end} = @(l) sqrt(max(1 - w(xi(l) + 1).^2, 0)) .* (xi(l) > -1);
            t  = nan(1, J + 2);
            Nf = J + 2;

        case 'mexhat'
            if tmax <= tmin
                error('graphfilterbank:range', ...
                    ['t_min = %.3g exceeds t_max = %.3g: the coarse end (sigma_max = %.1f mm) ' ...
                     'is finer than the fine end (sigma_min = %.1f mm).'], ...
                    tmin, tmax, 1000*sqrt(2*tmax), 1000*sqrt(2*tmin));
            end
            if isempty(Nf)
                nOct = log2(sqrt(tmax/tmin));
                Nf   = max(2, round(voices * nOct) + 1);
            end
            tw = logspace(log10(tmin), log10(tmax), Nf);
            g  = cell(1, Nf+1);
            % GSPBox's scaling function: quartic roll-off at 1.2x the wavelet peak, so it
            % fills the low-pass gap without dominating the max over (vertex, scale).
            g{1} = @(l) (1.2*exp(-1)) * exp(-(double(l(:))/(0.4*lminEff)).^4);
            for ii = 1:Nf, g{ii+1} = @(l) i_mexhat(l, tw(ii)); end
            t = [NaN, tw];

        case 'heat'
            % rheome.filters.frame pins the heat fine end at 1/lmax regardless of the tmin
            % policy, and runs one octave coarser at the top.
            tw = logspace(log10(2*tmax), log10(tmin), ...
                          i_heatcount(Nf, tmin, 2*tmax, voices));
            g = cell(1, numel(tw));
            for ii = 1:numel(tw), g{ii} = @(l) i_heat(l, tw(ii)); end
            t  = tw;
            Nf = numel(tw);

        otherwise
            error('graphfilterbank:wavelet', 'unknown family ''%s''.', family);
    end

    d = struct('g', {g}, 't', t, 'Nf', Nf);
end

function n = i_heatcount(Nf, tmin, tmax, voices)
    if ~isempty(Nf), n = Nf; return; end
    n = max(2, round(voices * log2(sqrt(tmax/tmin))) + 1);
end

function v = i_mexhat(l, t)
    x = double(l(:)) * t;
    v = x .* exp(-x);
end

function v = i_heat(l, t)
    v = exp(-(double(l(:)) * t));
end

% Author: Diellor Basha, 2026
