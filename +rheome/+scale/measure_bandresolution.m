function T = measure_bandresolution(name, S, snr)
% SCALE.MEASURE_BANDRESOLUTION  Resolution at each octave's OWN measured SNR, for one subject.
%
%   T = rheome.scale.measure_bandresolution(name)
%   T = rheome.scale.measure_bandresolution(name, rheome.scale.sensors(name), rheome.scale.bandsnr(name))
%
% The SNR per octave is this subject's (rheome.scale.bandsnr against its same-session noise run). SnrFixed = 10^(dB/20) (dB of power -> amplitude ratio);
% rheome.inverse.resolution on the left hemisphere, 200 seeds. Octaves with NaN SNR are skipped.
%
% Per band: snr_dB, SnrFixed, traceR, r50, lambda_loc, mtf_half, mtf_centroid, r50_superficial,
% r50_deep. Summary rows (band ""): r50_spread (max/min over bands), mtf_centroid_spread, and
% alpha_depth_ratio (deep/superficial r50 in 8-16 Hz).
% ⚠ For scale, a reference subject: r50 40.9-54.2 mm (x1.32), MTF centroid 140-157 mm (x1.12), alpha depth x2.59.
%
% See also: rheome.scale.bandsnr, rheome.inverse.resolution
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(S),   S = rheome.scale.sensors(name); end
    if nargin < 3 || isempty(snr), snr = rheome.scale.bandsnr(name); end
    SL = S.B.L.S;  T = rheome.scale.rows("bandresolution", strings(0,1), [], "");
    r50 = nan(height(snr),1);  cen = r50;  dr = r50;
    for i = 1:height(snr)
        dB = snr.snr_dB(i);  if isnan(dB), continue; end
        sf = 10^(dB/20);
        R = rheome.inverse.mne(S.G, S.ncm, struct('ChannelTypes',{S.chT}, 'InverseMeasure','amplitude', ...
                                           'nVert', S.nV, 'SnrFixed', sf));
        Lam = R.SNR / mean(R.SL.^2);  gs = (Lam*R.SL.^2) ./ (Lam*R.SL.^2 + 1);
        out = rheome.inverse.resolution(R.ImagingKernel, S.G, S.B.L.lbo, Vertices=SL.Vertices, Faces=SL.Faces, ...
                  Normals=SL.VertNormals, GlobalIdx=S.B.L.gv, SensorLoc=S.Loc, NumSeeds=200);
        q = discretize(out.depth, prctile(out.depth, [0 25 50 75 100]));
        md = accumarray(q, out.r50, [4 1], @median);
        r50(i) = 1e3*median(out.r50);  cen(i) = 1e3*out.mtfCentroid;  dr(i) = md(4)/md(1);
        T = [T; rheome.scale.rows("bandresolution", ...
            ["snr_dB" "SnrFixed" "traceR" "r50" "lambda_loc" "mtf_half" "mtf_centroid" "r50_superficial" "r50_deep"], ...
            [dB sf sum(gs) r50(i) 1e3*out.lambdaLoc 1e3*out.mtfHalf cen(i) 1e3*md(1) 1e3*md(4)], ...
            ["dB" "ratio" "modes" "mm" "mm" "mm" "mm" "mm" "mm"], snr.band(i))]; %#ok<AGROW>
    end
    ia = find(snr.band == "8-16 Hz", 1);
    T = [T; rheome.scale.rows("bandresolution", ["r50_spread" "mtf_centroid_spread" "alpha_depth_ratio"], ...
        [max(r50)/min(r50) max(cen)/min(cen) dr(ia)], "ratio")];
end

% Author: Diellor Basha, 2026
