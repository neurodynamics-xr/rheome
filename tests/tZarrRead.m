classdef tZarrRead < matlab.unittest.TestCase
% rheome.io.read.zarr against arrays written by zarr-python (tests/fixtures/zarr/make.py; values arange in C order): one chunk,
% a grid with ragged edge chunks, 1-D, shards (an absent inner chunk read as the fill value, the index at the end or the start),
% 3-D shards, partial ranges across chunk and shard edges; the codecs nxr writes (blosc whole and as the inner chunks of shards) and
% zstd and gzip, decoded by numcodecs through MATLAB's Python ($RHEOME_PYTHON names a python (3.9-3.11, numcodecs installed) for
% pyenv, if MATLAB's Python has none; those tests are skipped when numcodecs cannot be imported), and a transposed array refused.
%
% Author: Diellor Basha, 2026

    properties (TestParameter)
        name = {'single', 'grid', 'vector', 'shard', 'shard3', 'shardstart'}
        packed = {'compressed', 'blosc', 'blosc_shard', 'gzip', 'zstd'}
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

        function compressedIsDecoded(tc, packed)
            i_numcodecs(tc);
            [x, meta] = rheome.io.read.zarr(i_dir(packed));
            tc.verifyEqual(x, i_expect(packed, meta));
        end

        function compressedRangesAcrossChunkEdges(tc, packed)
            i_numcodecs(tc);
            [full, meta] = rheome.io.read.zarr(i_dir(packed));
            shape = meta.shape(:)'; rng = RandStream('mt19937ar', 'Seed', 2);
            for t = 1:10
                r = cell(1, numel(shape));
                for d = 1:numel(shape), r{d} = sort(randi(rng, shape(d), 1, 2)); end
                idx = cellfun(@(a) a(1):a(2), r, 'UniformOutput', false);
                tc.verifyEqual(rheome.io.read.zarr(i_dir(packed), r), full(idx{:}), sprintf('%s range %s', packed, mat2str(cell2mat(r))));
            end
        end

        function missingNumcodecsSaysWhatToInstall(tc)
            % without numcodecs the error names it, and pyenv
            tc.assumeFalse(i_has_numcodecs(), 'numcodecs is importable here');
            tc.verifyError(@() rheome.io.read.zarr(i_dir('blosc')), 'io:read:zarr:codec');
        end

        function transposedIsRefused(tc)
            tc.verifyError(@() rheome.io.read.zarr(i_dir('transposed')), 'io:read:zarr:layout');
        end

        function rangeOutsideIsRefused(tc)
            tc.verifyError(@() rheome.io.read.zarr(i_dir('grid'), {[], [4 6]}), 'io:read:zarr:range');
        end
    end
end

function i_numcodecs(tc)
tc.assumeTrue(i_has_numcodecs(), 'numcodecs is not importable in MATLAB''s Python (set RHEOME_PYTHON)');
end

function ok = i_has_numcodecs()
v = getenv('RHEOME_PYTHON');
if ~isempty(v) && pyenv().Status == "NotLoaded", pyenv(Version=v); end
try, py.importlib.import_module('numcodecs'); ok = true; catch, ok = false; end
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
if ismember(name, {'shard', 'blosc_shard'}), x(1, 17:32) = -1; end
end
