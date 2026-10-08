function Wcv = curlvec(ctx)
% FLOW.CURLVEC  Ambient curl-VECTOR kernel [3V x C]  (the full (grad x J) per vertex).
%   Wcv = rheome.flow.curlvec(ctx)
%
% Unlike rheome.flow.curl (which projects the curl onto the surface normal to give scalar vorticity),
% this keeps the full ambient 3-vector curl, area-weighted to vertices, interleaved
% [x1 y1 z1 x2 ...]. Used by rheome.flow.helicity (J . (grad x J)).
%
% See also: rheome.flow.curl, rheome.flow.helicity
%
% Author: Diellor Basha, 2026

    fg = ctx.fg;  currentKernel = ctx.currentKernel;
    Jx = currentKernel(1:3:end,:); Jy = currentKernel(2:3:end,:); Jz = currentKernel(3:3:end,:);
    cvx = fg.Gy*Jz - fg.Gz*Jy;      % (grad x J)_x per face  [nF x C]
    cvy = fg.Gz*Jx - fg.Gx*Jz;
    cvz = fg.Gx*Jy - fg.Gy*Jx;
    Vx = fg.W*cvx; Vy = fg.W*cvy; Vz = fg.W*cvz;          % area-weighted to vertices [V x C]
    Wcv = zeros(3*fg.nV, size(currentKernel,2));
    Wcv(1:3:end,:) = Vx; Wcv(2:3:end,:) = Vy; Wcv(3:3:end,:) = Vz;
end

% Author: Diellor Basha, 2026
