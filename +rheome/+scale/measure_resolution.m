function T = measure_resolution(name, S)
% SCALE.MEASURE_RESOLUTION  The resolution report's headline numbers for one subject.
%
%   T = rheome.scale.measure_resolution(name)
%   T = rheome.scale.measure_resolution(name, rheome.scale.sensors(name))
%
% Same call as resolution_scales_omega.m:52: rheome.inverse.resolution on the left hemisphere, the plain
% MNE at SnrFixed = 3, 300 seeds. Returns rheome.scale.rows of
%   r50_median, r50_q25, r50_q75   geodesic radius holding half the PSF power        mm
%   ple_median                      peak localisation error                           mm
%   offhemi_median                  fraction of PSF power on the other hemisphere     fraction
%   mtf_half, mtf_centroid          MTF half-max (⚠ jumpy) and centroid (⭐)           mm
%   lambda_loc                      the localisation floor as a wavelength            mm
%   dof_subspace                    sum of the MTF                                     modes
%   r50_depth_q1 .. _q4             median r50 by depth quartile (superficial->deep)   mm
%   depth_rho                       Spearman rho(depth, r50)
%
% ⚠ For scale: a reference subject (270 CTF channels) reads r50 52 mm, PLE 23 mm, 12 % off-hemisphere,
% MTF half-max 99 mm.
%
% See also: rheome.inverse.resolution, rheome.scale.sensors
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(S), S = rheome.scale.sensors(name); end
    SL = S.B.L.S;
    out = rheome.inverse.resolution(S.Res.ImagingKernel, S.G, S.B.L.lbo, Vertices=SL.Vertices, Faces=SL.Faces, ...
              Normals=SL.VertNormals, GlobalIdx=S.B.L.gv, SensorLoc=S.Loc, NumSeeds=300);
    q = discretize(out.depth, prctile(out.depth, [0 25 50 75 100]));
    md = accumarray(q, out.r50, [4 1], @median);
    m  = ["r50_median" "r50_q25" "r50_q75" "ple_median" "offhemi_median" "mtf_half" ...
          "mtf_centroid" "lambda_loc" "dof_subspace" "r50_depth_q1" "r50_depth_q2" ...
          "r50_depth_q3" "r50_depth_q4" "depth_rho"];
    v  = [1e3*median(out.r50) 1e3*prctile(out.r50,25) 1e3*prctile(out.r50,75) 1e3*median(out.ple) ...
          median(out.offHemi) 1e3*out.mtfHalf 1e3*out.mtfCentroid 1e3*out.lambdaLoc out.dofSubspace ...
          1e3*md(:)' corr(out.depth(:), out.r50(:), 'Type', 'Spearman')];
    u  = ["mm" "mm" "mm" "mm" "fraction" "mm" "mm" "mm" "modes" "mm" "mm" "mm" "mm" "rho"];
    T = rheome.scale.rows("resolution", m, v, u);
end

% Author: Diellor Basha, 2026
