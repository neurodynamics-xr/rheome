function [U, defined] = sen_pattern(arr, pl, name, lam, f0, tv)
% ⭐ ONE TABLE OF GENERATORS, shared by the family figure and the web export, so the two can
% never drift apart. Every pattern here is PERIODIC at f0, which is what lets an animation
% loop seamlessly; diffusion is not, and appears only in the static figure.
%
% ⚠ ROTATION, SPIRAL AND TARGET ARE UNDEFINED ON A CHAIN. They are built from the polar
% angle about the array centre, and a 1D array has no such angle -- every sensor sits at
% phi = 0 or pi. `defined` says so rather than returning a degenerate field that looks like
% a measurement.
    k = 2*pi/lam;  w = 2*pi*f0;
    kh = sen_dir(pl);  x = pl.Pc * kh(:);
    phi = atan2(pl.x2, pl.x1);
    r   = hypot(pl.x1, pl.x2);
    % ⚠ THE RADIAL WAVENUMBER IS THE REQUESTED ONE, not one read off the array. Deriving it
    % from max(r) made spiral and target IGNORE lam entirely, so a scale control moved the
    % travelling wave and left the rings untouched -- a control that silently does nothing
    % on half the generators.
    kr  = k;
    defined = true;
    needs2D = any(strcmp(name, {'rotating','spiral','target'}));
    if needs2D && arr.Dim < 2
        defined = false;  U = zeros(arr.nCh, numel(tv));  return;
    end
    switch name
        case 'travelling', f = @(P,t) cos(k*x - w*t);
        case 'standing',   f = @(P,t) cos(k*x) .* cos(w*t);
        % ⚠ A SINGLE-ARM ROTATION HAS NO WAVELENGTH. cos(phi - w*t) is one crest sweeping
        % once around the centre, and its spatial content is set by the array's own extent,
        % not by lam. It is the one generator here that a scale control cannot move, and
        % saying so beats pretending otherwise.
        case 'rotating',   f = @(P,t) cos(phi - w*t);
        case 'spiral',     f = @(P,t) cos(phi + kr*r - w*t);
        case 'target',     f = @(P,t) cos(kr*r - w*t);
        case 'colliding'
            % Two fronts crossing at 90 degrees: a moving interference lattice, which is
            % neither a clean travelling wave nor a standing one.
            kh2 = ([-kh(2) kh(1) 0]);  x2 = pl.Pc * kh2(:);
            f = @(P,t) 0.5*(cos(k*x - w*t) + cos(k*x2 - w*t));
        otherwise, error('demos:pattern:name', 'unknown pattern ''%s''.', name);
    end
    U = rheome.sensors.sample(arr, f, tv);
end

% Author: Diellor Basha, 2026
