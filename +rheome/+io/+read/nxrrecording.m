function R = nxrrecording(store, subject, recording, opts)
% IO.READ.NXRRECORDING  Read a sensor recording, or one span of its samples, from an nxr store.
%
%   R = rheome.io.read.nxrrecording(store, subject, recording)                       % metadata, channels and flags; no samples
%   R = rheome.io.read.nxrrecording(store, subject, recording, Samples=[1 32768])    % and those samples of every channel
%
% A recording is a field of kind 'recording' on the product (channels x a time Line): its channels a points manifold whose
% members are named by the partition with one code a member (<set>_names), its time Line a lattice whose dimension row states
% the spacing (1/fs) and origin (t0) under the 'clock' coordinate system. Its bad channels are the <recording>_flags labelling
% ('bad'). The lookups are the editor's (adapters/nxr recordingFromRows, lineFromRows, setFromRows).
%
% INPUTS:  store the dataset folder (holds dataset.sqlite);  subject its subject.name;  recording the field's name
% OPTIONS: Samples [first last] 1-based, inclusive (default: none read)
%
% OUTPUT (struct R)
%   .Data      [C x n] single, the samples asked for ([] when none), rows in channel order, as stored
%   .Fs        sampling rate (Hz) = 1/spacing;  .T0  the time of sample 1 (s);  .nSamples  the whole recording's length
%   .Channels  [C x 1] string, the members' names;  .Bad [C x 1] logical (all false when no flags are stored)
%   .Unit      the field's unit ('' when not stated);  .Samples the span read;  .Name, .Session, .Id, .Path
%
% See also: rheome.io.read.nxrkernel, rheome.io.read.zarr
%
% Author: Diellor Basha, 2026
arguments, store (1,1) string, subject (1,1) string, recording (1,1) string, opts.Samples double = [], end
q = @(s) "'" + replace(s, "'", "''") + "'";
f = nxrsql(store, "SELECT f.id, f.name, f.path, f.session, f.unit, f.manifold_id, u.path AS store FROM field f JOIN subject u ON u.id = f.subject_id " + ...
    "WHERE u.name = " + q(subject) + " AND u.status = 'complete' AND f.status = 'complete' AND f.kind = 'recording' AND f.name = " + q(recording));
if isempty(f), error('io:read:nxr:notFound', '%s: no recording ''%s''', subject, recording); end
fac = nxrsql(store, "SELECT m.id, m.type, m.n_vertices, m.geometry_id FROM geometry_dimension d JOIN manifold p ON p.geometry_id = d.geometry_id " + ...
    "JOIN manifold m ON m.id = d.manifold_id WHERE p.id = " + q(f.manifold_id) + " ORDER BY d.ordinal");
if numel(fac) < 2 || ~strcmp(fac(2).type, 'line'), error('io:read:nxr:layout', 'recording ''%s'': not on channels x a time Line', recording); end
ln = nxrsql(store, "SELECT d.spacing, d.origin, c.name AS cs FROM geometry_dimension d JOIN geometry g ON g.id = d.geometry_id " + ...
    "LEFT JOIN coordinate_system c ON c.id = g.coordinate_system_id WHERE d.geometry_id = " + q(fac(2).geometry_id) + " AND d.ordinal = 0");
if isempty(ln) || isempty(ln.spacing) || ~strcmp(ln.cs, 'clock'), error('io:read:nxr:layout', 'recording ''%s'': its time Line states no clock spacing', recording); end
z = fullfile(store, f.store);
C = fac(1).n_vertices; T = fac(2).n_vertices;
R = struct('Data', single([]), 'Fs', 1 / ln.spacing, 'T0', i_num(ln.origin), 'nSamples', T, 'Channels', i_names(store, z, fac(1).id, C, q), ...
    'Bad', false(C, 1), 'Unit', i_str(f.unit), 'Samples', [], 'Name', f.name, 'Session', i_str(f.session), 'Id', f.id, 'Path', fullfile(z, f.path));
fl = nxrsql(store, "SELECT s.array_path, s.dictionary_id FROM selection s WHERE s.manifold_id = " + q(fac(1).id) + " AND s.name = " + q(f.name + "_flags") + ...
    " AND s.array_path IS NOT NULL AND s.dictionary_id IS NOT NULL AND s.status = 'complete'");
if ~isempty(fl)
    bad = nxrsql(store, "SELECT code FROM dictionary_entry WHERE dictionary_id = " + q(fl(1).dictionary_id) + " AND name = 'bad'");
    if ~isempty(bad), R.Bad = ismember(double(rheome.io.read.zarr(fullfile(z, fl(1).array_path))), [bad.code]); end
end
if ~isempty(opts.Samples)
    R.Samples = opts.Samples([1 end]);
    R.Data = rheome.io.read.zarr(R.Path, {[], R.Samples});
end
end

function names = i_names(store, z, mid, n, q)
% the members named by the partition with one code a member: <set>_names first, else one with a level a member (setFromRows)
s = nxrsql(store, "SELECT name, array_path, dictionary_id, n_members FROM selection WHERE manifold_id = " + q(mid) + " AND type = 'set' " + ...
    "AND array_path IS NOT NULL AND dictionary_id IS NOT NULL AND member_of_id IS NULL AND n_elements = " + n + " AND status = 'complete'");
names = "channel " + (1:n)';
if isempty(s), return, end
k = find(endsWith(string({s.name}), "_names"), 1);
if isempty(k), k = find([s.n_members] == n, 1); end
if isempty(k), return, end
codes = double(rheome.io.read.zarr(fullfile(z, s(k).array_path)));
e = nxrsql(store, "SELECT code, name FROM dictionary_entry WHERE dictionary_id = " + q(s(k).dictionary_id));
[ok, at] = ismember(codes, [e.code]);
names(ok) = string({e(at(ok)).name});
end

function v = i_num(x)
if isempty(x), v = 0; else, v = x; end
end

function v = i_str(x)
if isempty(x), v = ''; else, v = char(x); end
end
