% Author: Diellor Basha, 2026
function tf = isemptyroom(name)
% ISEMPTYROOM  True for an empty-room (noise) MEG recording, from its BIDS name.
%
%   tf = rheome.qc.isemptyroom(name)
%
%   name : char/string, or a cellstr/string array of them -- a recording name, run stem or path
%          (e.g. 'sub-01_ses-04_task-noise_run-02_meg', a data_0raw_*.mat path, a *_meg.ds path).
%   tf   : logical, the same size as name.
%
%   An empty-room recording is one whose BIDS entities say so: 'task-noise' (used by some open MEG datasets),
%   'task-emptyroom' (the BIDS convention), or subject 'sub-emptyroom'. Matching is on whole
%   entities (bounded by '_', '/', '\' or the end), case-insensitive, so 'task-noisetest' or a
%   folder that merely contains the word 'noise' does not count.
%   ⚠ Detector v2 (PROTOCOL-v2): an empty-room recording holds only sensor noise, so its spectrum
%   is flat by construction and no bad-channel rule may label it. Decided by Diellor Basha, 2026-10-02.
%
% See also: rheome.qc.badchannels

pat = '(^|[_/\\])(task-(noise|emptyroom)|sub-emptyroom)(?=[_/\\.]|$)';
if ischar(name), tf = ~isempty(regexpi(name, pat, 'once')); return; end
name = cellstr(name);
tf = reshape(~cellfun(@isempty, regexpi(name, pat, 'once')), size(name));
end
% Author: Diellor Basha, 2026
