function [L, T, M] = sel_labelio(db, L, T, M)
% SEL_LABELIO  Read (nargin 1) or atomically write (nargin 3 or 4) the note tables.
%   Written to a temporary name beside the file and renamed: a reader never sees a
%   partial write, and a crash leaves the previous state intact -- the transaction.
%
% Three tables share one sidecar file and one transaction log: `label` (categorical marks),
% `measure` (numbers measured at a node) and `txn`. Writing fewer than three keeps the
% others as they are, so a label transaction cannot drop measurements.
%
% ⚠ A FILE WRITTEN BEFORE `measure` EXISTED HAS ONLY TWO TABLES. Loading tolerates that and
% returns an empty measure table rather than erroring, so old stores keep working.
%
% Author: Diellor Basha, 2026

    f = db.labelFile;
    if nargin == 1
        if exist(f, 'file') == 2
            s = load(f, 'label', 'txn');  L = s.label;  T = s.txn;
            w = whos('-file', f);
            if ismember('measure', {w.name}), sm = load(f, 'measure');  M = sm.measure; else, M = i_empty(); end
        else
            M = i_empty();
            L = table('Size', [0 10], 'VariableTypes', {'double','string','double','double','double','double','string','double','string','double'}, ...
                      'VariableNames', {'label_id','recording_id','channel_id','level','k','band_id','kind','value','source','txn_id'});
            T = table('Size', [0 7], 'VariableTypes', {'double','string','string','double','string','string','string'}, ...
                      'VariableNames', {'txn_id','recording_id','kind','count','content_hash','author','created'});
        end
        return
    end
    if nargin < 4, [~, ~, M] = sel_labelio(db); end
    label = L;  txn = T;  measure = M;                          %#ok<NASGU>
    tmp = [f '.tmp'];
    save(tmp, 'label', 'txn', 'measure', '-v7');
    movefile(tmp, f, 'f');
end

function M = i_empty()
    M = table('Size', [0 11], ...
              'VariableTypes', {'double','string','string','double','double','double','double','string','double','string','double'}, ...
              'VariableNames', {'measure_id','recording_id','scope','unit_id','level','k','band_id','kind','value','source','txn_id'});
end
% Author: Diellor Basha, 2026
