classdef tZarrRead < matlab.unittest.TestCase
% rheome.io.read.zarr against arrays written by zarr-python (tests/fixtures/zarr/make.py; values arange in C order): one chunk,
% a grid with ragged edge chunks, 1-D, shards (an absent inner chunk read as the fill value, the index at the end or the start),
% 3-D shards, partial ranges across chunk and shard edges, and a compressed array refused.
%
% Author: Diellor Basha, 2026

    properties (TestParameter)
        name = {'single', 'grid', 'vector', 'shard', 'shard3', 'shardstart'}
    end

    methods (Test)
        function wholeArray(tc, name)
            [x, meta] = rheome.io.read.zarr(i_dir(name));
            tc.verifyEqual(x, i_expect(name, meta));
        end

        function everyRangeAcrossChunkEdges(tc, name)
            [full, meta] = rheome.io.read.zarr(i_dir(name));
            shape = meta.shape(:)'; rng = RandStream('mt19937ar', 'Seed', 1);
            for t = 1:25
                r = cell(1, numel(shape));
                for d = 1:numel(shape), r{d} = sort(randi(rng, shape(d), 1, 2)); end
                idx = cellfun(@(a) a(1):a(2), r, 'UniformOutput', false);
                if isscalar(shape), idx{2} = 1; end
                tc.verifyEqual(rheome.io.read.zarr(i_dir(name), r), full(idx{:}), sprintf('%s range %s', name, mat2str(cell2mat(r))));
            end
        end

        function absentInnerChunkIsTheFill(tc)
            x = rheome.io.read.zarr(i_dir('shard'), {[1 1], [17 32]});
            tc.verifyEqual(x, -ones(1, 16, 'single'));
        end

        function compressedIsRefused(tc)
            tc.verifyError(@() rheome.io.read.zarr(i_dir('compressed')), 'io:read:zarr:layout');
        end

        function rangeOutsideIsRefused(tc)
            tc.verifyError(@() rheome.io.read.zarr(i_dir('grid'), {[], [4 6]}), 'io:read:zarr:range');
        end
    end
end

function d = i_dir(name)
d = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'zarr', name);
end

function x = i_expect(name, meta)
% arange in C order, of the array's class; the 'shard' fixture's hole (row 0, columns 16..31) holds its fill value -1
shape = meta.shape(:)'; nd = numel(shape);
cls = struct('float64', 'double', 'float32', 'single', 'int32', 'int32', 'uint8', 'uint8').(meta.data_type);
x = cast(0:prod(shape) - 1, cls);
if nd > 1, x = permute(reshape(x, fliplr(shape)), nd:-1:1); else, x = x(:); end
if strcmp(name, 'shard'), x(1, 17:32) = -1; end
end
