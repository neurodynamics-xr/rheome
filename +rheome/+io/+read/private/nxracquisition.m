function a = nxracquisition(store, id)
% NXRACQUISITION  An acquisition row of an nxr store (schema 50, D161) as a struct: id, label and the BIDS entities ses, task, run,
% acq, proc, split (each '' where the source does not say). Empty id: the same fields, all ''.
a = struct('id', '', 'label', '', 'ses', '', 'task', '', 'run', '', 'acq', '', 'proc', '', 'split', '');
if isempty(id), return, end
r = nxrsql(store, "SELECT id, label, ses, task, run, acq, proc, split FROM acquisition WHERE id = '" + replace(string(id), "'", "''") + "'");
if isempty(r), return, end
for f = string(fieldnames(a))'
    if ~isempty(r.(f)), a.(f) = char(r.(f)); end
end
end
