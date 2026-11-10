function [x, meta] = zarr(folder, range)
% IO.READ.ZARR  Read an uncompressed zarr v3 array, or any rectangular part of it, without the zarr toolbox.
%
%   x = rheome.io.read.zarr(folder)                       % the whole array
%   x = rheome.io.read.zarr(folder, {[], [s0 s1]})        % every row, columns s0..s1 (1-based, inclusive)
%   [x, meta] = rheome.io.read.zarr(folder, {[r0 r1], []})
%
% folder   the array's folder (holds its zarr.json)
% range    one [first last] per dimension, 1-based and inclusive; [] for the whole dimension
% x        the part, in the array's own class (float32 -> single), indexed as the zarr array is: x(i,j) = a[i-1, j-1]
%          (a 1-D array is a column)
% meta     the zarr.json (shape, data_type, chunk_grid, codecs, attributes, ...)
%
% Reads three layouts: one chunk, a regular grid of chunks, and sharding_indexed shards (each shard read only for the inner
% chunks the range touches, by its index). A chunk or inner chunk absent from disk is the fill value. Only raw little-endian
% bytes are read: a compressed or transposed array is refused, never misread.
% shortcut: the shard index's crc32c is skipped, not checked; check it if shards ever arrive from outside the hub's writers.
%
% Example:
%   x = rheome.io.read.zarr(fullfile(store, rec), {[], [1 32768]});      % one 2^15-sample chunk of every channel
%
% See also: rheome.io.read.nxrsurface
%
% Author: Diellor Basha, 2026
arguments, folder (1,1) string, range cell = {}, end
meta = jsondecode(fileread(fullfile(folder, 'zarr.json')));
if ~isfield(meta, 'zarr_format') || meta.zarr_format ~= 3 || ~strcmp(meta.node_type, 'array')
    error('io:read:zarr:layout', '%s: not a zarr v3 array', folder);
end
types = struct('float64', 'double', 'float32', 'single', 'int8', 'int8', 'int16', 'int16', 'int32', 'int32', 'int64', 'int64', ...
    'uint8', 'uint8', 'uint16', 'uint16', 'uint32', 'uint32', 'uint64', 'uint64', 'bool', 'uint8');
if ~isfield(types, meta.data_type), error('io:read:zarr:layout', '%s: data type %s is not read', folder, meta.data_type); end
cls = types.(meta.data_type);
shape = meta.shape(:)'; nd = numel(shape);
outer = meta.chunk_grid.configuration.chunk_shape(:)';
if ~strcmp(meta.chunk_grid.name, 'regular'), error('io:read:zarr:layout', '%s: chunk grid %s', folder, meta.chunk_grid.name); end
codecs = i_list(meta.codecs);
inner = [];
if numel(codecs) == 1 && strcmp(codecs{1}.name, 'sharding_indexed')
    sc = codecs{1}.configuration;
    i_raw(i_list(sc.codecs), folder);
    ic = i_list(sc.index_codecs);
    i_raw(ic(1), folder);
    if numel(ic) > 2 || (numel(ic) == 2 && ~strcmp(ic{2}.name, 'crc32c')), error('io:read:zarr:layout', '%s: shard index codecs', folder); end
    inner = sc.chunk_shape(:)'; idxCrc = 4 * (numel(ic) == 2); atEnd = ~isfield(sc, 'index_location') || strcmp(sc.index_location, 'end');
else
    i_raw(codecs, folder);
end
if isempty(range), range = cell(1, nd); end
if numel(range) ~= nd, error('io:read:zarr:range', '%s: %d ranges for %d dimensions', folder, numel(range), nd); end
lo = ones(1, nd); hi = shape;
for d = 1:nd
    if ~isempty(range{d}), lo(d) = range{d}(1); hi(d) = range{d}(end); end
    if lo(d) < 1 || hi(d) > shape(d) || lo(d) > hi(d), error('io:read:zarr:range', '%s: range %d..%d outside 1..%d', folder, lo(d), hi(d), shape(d)); end
end
fill = cast(i_fill(meta.fill_value), cls);
x = repmat(fill, [hi - lo + 1, 1]);                              % [n 1] for a 1-D array: a column
key = @(c) i_key(meta.chunk_key_encoding, c);
for c = i_cells(lo, hi, outer)'                                  % the outer chunks (shards) the range touches
    c = c'; o = c .* outer;                                      %#ok<FXSET> 0-based origin of the chunk
    f = fullfile(folder, key(c));
    if isempty(inner)
        if ~isfile(f), continue, end
        x = i_place(x, i_bytes(f, 0, prod(outer), cls, folder), outer, o, lo, hi);
        continue
    end
    if ~isfile(f), continue, end
    nin = outer ./ inner;
    fid = fopen(f, 'r', 'ieee-le'); cl = onCleanup(@() fclose(fid));
    nIdx = prod(nin) * 16;
    fseek(fid, 0, 'eof'); len = ftell(fid);
    if atEnd, fseek(fid, len - nIdx - idxCrc, 'bof'); else, fseek(fid, 0, 'bof'); end
    ix = fread(fid, [2, prod(nin)], 'uint64=>uint64');          % C order over the inner grid: (offset, nbytes) per inner chunk
    ilo = max(lo, o + 1); ihi = min(hi, o + outer);              % the part of the range in this shard
    stride = fliplr(cumprod(fliplr([nin(2:end) 1])));            % C order over the inner grid
    for q = i_cells(ilo - o, ihi - o, inner)'
        q = q'; k = sum(q .* stride) + 1;                        %#ok<FXSET>
        if ix(1,k) == intmax('uint64'), continue, end            % absent: the fill value
        if double(ix(2,k)) ~= prod(inner) * i_size(cls), error('io:read:zarr:layout', '%s: inner chunk of %d bytes, not raw', f, ix(2,k)); end
        x = i_place(x, i_bytes(fid, double(ix(1,k)), prod(inner), cls, f), inner, o + q .* inner, lo, hi);
    end
    clear cl
end
if strcmp(meta.data_type, 'bool'), x = x ~= 0; end
end

function x = i_place(x, raw, cs, o, lo, hi)
% a chunk's values (C order, chunk shape cs, 0-based origin o) into the part x covering lo..hi
nd = numel(cs);
if nd > 1, b = permute(reshape(raw, fliplr(cs)), nd:-1:1); else, b = raw(:); end
src = cell(1, max(nd, 2)); dst = src; src(:) = {1}; dst(:) = {1};
for d = 1:nd
    g = max(lo(d), o(d) + 1):min(hi(d), o(d) + cs(d));          % global 1-based indices covered
    src{d} = g - o(d); dst{d} = g - lo(d) + 1;
end
x(dst{:}) = b(src{:});
end

function v = i_bytes(f, offset, n, cls, name)
if isnumeric(f), fid = f; else, fid = fopen(f, 'r', 'ieee-le'); cl = onCleanup(@() fclose(fid)); end
fseek(fid, offset, 'bof');
v = fread(fid, n, [cls '=>' cls]);
if numel(v) ~= n, error('io:read:zarr:layout', '%s: %d of %d values (a compressed or truncated chunk)', name, numel(v), n); end
end

function cells = i_cells(lo, hi, cs)
% the 0-based chunk coordinates (one row each) covering lo..hi (1-based)
g = arrayfun(@(d) floor((lo(d) - 1) / cs(d)):floor((hi(d) - 1) / cs(d)), 1:numel(cs), 'UniformOutput', false);
[g{:}] = ndgrid(g{:});
cells = cell2mat(cellfun(@(a) a(:), g, 'UniformOutput', false));
end

function k = i_key(enc, c)
sep = '/'; if isfield(enc, 'configuration') && isfield(enc.configuration, 'separator'), sep = enc.configuration.separator; end
parts = strjoin(string(c), sep);
switch enc.name
    case 'default', k = "c" + sep + parts;
    case 'v2', k = parts;
    otherwise, error('io:read:zarr:layout', 'chunk key encoding %s', enc.name);
end
end

function i_raw(codecs, folder)
% the only array codec read: raw bytes, little-endian
if numel(codecs) ~= 1 || ~strcmp(codecs{1}.name, 'bytes')
    error('io:read:zarr:layout', '%s: codecs %s are not read (compressed or transposed): refused rather than misread', folder, strjoin(cellfun(@(c) c.name, codecs, 'UniformOutput', false), '+'));
end
if isfield(codecs{1}, 'configuration') && isfield(codecs{1}.configuration, 'endian') && ~strcmp(codecs{1}.configuration.endian, 'little')
    error('io:read:zarr:layout', '%s: big-endian', folder);
end
end

function c = i_list(c)
% jsondecode gives a struct array for codecs of one shape, a cell for mixed: always a cell
if isstruct(c), c = num2cell(c); end
end

function v = i_fill(f)
if isnumeric(f) || islogical(f), v = double(f); return, end
switch string(f), case "NaN", v = NaN; case "Infinity", v = Inf; case "-Infinity", v = -Inf; otherwise, v = 0; end
end

function n = i_size(cls)
n = numel(typecast(zeros(1, 1, cls), 'uint8'));
end
