% Author: Diellor Basha, 2026
function [flag, why, rawFlag, rawWhy] = badchannelrules(M, P)
% BADCHANNELRULES  The bad-channel rules of rheome.qc.badchannels and how they combine into a label (protocol v3).
%
%   [flag, why, rawFlag, rawWhy] = rheome.qc.badchannelrules(M, P)
%
%   M : metrics from rheome.qc.badchannels (fields sd, z_sd, r_sd, jumpfrac, z_lo, r_lo, z_hi, r_hi), each [nC x 1].
%   P : thresholds (rheome.qc.badchannels i_params: zFlat rFlat zNoisy rNoisy jumpFrac zPsd lrPsd).
%   flag    [nC x 1] logical : the channel is LABELLED bad.
%   why     {nC x 1}         : every rule that fired on a labelled channel ('' otherwise), e.g. 'NOISY PSD_HIGH'.
%   rawFlag [nC x 1] logical : ANY rule fired, labelled or not -- kept for the record.
%   rawWhy  {nC x 1}         : every rule that fired, whether or not the channel is labelled.
%
%   ⭐ v3 (2026-10-08): PSD_HIGH (65–115 Hz) never labels a channel on its own; it
%   only corroborates FLAT, NOISY, JUMPY or PSD_LOW. ⚠ Across one resting cohort all 53 PSD_HIGH-only subject marks of v1/v2
%   were sensor noise-floor differences above 65 Hz (38 on MRO11/MRO12/MRP41/MRP51/MRP53, each also raised in
%   the same session's empty-room run), not bad signal in the band we analyse.
%   The empty-room gate (rheome.qc.isemptyroom) is applied by the caller, after this step.
%
% See also: rheome.qc.badchannels, rheome.qc.isemptyroom

n = numel(M.sd);
R.FLAT  = M.sd == 0 | (M.z_sd < P.zFlat & M.r_sd < P.rFlat);
R.NOISY = M.z_sd > P.zNoisy & M.r_sd > P.rNoisy;
R.JUMPY = M.jumpfrac > P.jumpFrac;
R.PSD_LOW  = abs(M.z_lo) > P.zPsd & abs(log10(M.r_lo)) > P.lrPsd;
R.PSD_HIGH = abs(M.z_hi) > P.zPsd & abs(log10(M.r_hi)) > P.lrPsd;
fn = fieldnames(R); rawFlag = false(n, 1); rawWhy = repmat({''}, n, 1);
for i = 1:numel(fn)
    hit = R.(fn{i}); rawFlag = rawFlag | hit(:);
    for c = find(hit(:))', rawWhy{c} = strtrim([rawWhy{c} ' ' fn{i}]); end
end
flag = rawFlag & ~(R.PSD_HIGH(:) & ~(R.FLAT(:) | R.NOISY(:) | R.JUMPY(:) | R.PSD_LOW(:)));
why = rawWhy; why(~flag) = {''};
end
% Author: Diellor Basha, 2026
