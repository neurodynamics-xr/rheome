"""Writes the zarr fixtures of tZarrRead with zarr-python (the reference writer): values arange in C order, synthetic only.

    uv run --frozen --project <nxr checkout> python -I make.py      (from this folder; rewrites every array)
"""
import shutil
from pathlib import Path

import numpy as np
import zarr

here = Path(__file__).parent


def put(name, shape, dtype, chunks, shards=None, fill=0, compressors=None, index_location="end", hole=None):
    p = here / name
    shutil.rmtree(p, ignore_errors=True)
    a = np.arange(np.prod(shape)).reshape(shape).astype(dtype)
    if hole is not None:
        a[hole] = fill                                   # equal to the fill value: zarr leaves that inner chunk off disk
    kw = {"shards": shards} if shards else {}
    if shards and index_location != "end":
        from zarr.codecs import BytesCodec, ShardingCodec
        z = zarr.create_array(store=str(p), shape=shape, dtype=dtype, chunks=shards, fill_value=fill, zarr_format=3, compressors=None,
                              serializer=ShardingCodec(chunk_shape=chunks, codecs=[BytesCodec()], index_location=index_location), filters=None)
    else:
        z = zarr.create_array(store=str(p), shape=shape, dtype=dtype, chunks=chunks, fill_value=fill, zarr_format=3, compressors=compressors, **kw)
    z[...] = a


put("single", (5, 3), "float64", (5, 3))
put("grid", (7, 5), "int32", (3, 2))
put("vector", (11,), "uint8", (4,))
put("shard", (6, 100), "float32", (1, 16), shards=(4, 64), fill=-1, hole=(0, slice(16, 32)))
put("shard3", (3, 4, 10), "float64", (1, 2, 4), shards=(2, 4, 8))
put("shardstart", (6, 100), "float32", (1, 16), shards=(4, 64), index_location="start")
put("compressed", (4, 4), "float64", (4, 4), compressors="auto")
# the codecs nxr writes (blosc: zstd, level 3, byte shuffle; arrays over 4 MiB), whole and as the inner chunks of shards, and the others zarr names
from zarr.codecs import BloscCodec, GzipCodec, TransposeCodec, ZstdCodec
blosc = BloscCodec(cname="zstd", clevel=3, shuffle="shuffle")
put("blosc", (7, 5), "float32", (3, 2), compressors=blosc)
put("blosc_shard", (6, 100), "float32", (1, 16), shards=(4, 64), fill=-1, hole=(0, slice(16, 32)), compressors=blosc)
put("gzip", (7, 5), "int32", (3, 2), compressors=GzipCodec(level=5))
put("zstd", (7, 5), "float64", (3, 2), compressors=ZstdCodec(level=3))
p = here / "transposed"                                  # a codec the reader does not do: refused
shutil.rmtree(p, ignore_errors=True)
z = zarr.create_array(store=str(p), shape=(4, 3), dtype="float64", chunks=(4, 3), zarr_format=3, filters=[TransposeCodec(order=(1, 0))], compressors=None)
z[...] = np.arange(12.0).reshape(4, 3)
print(sorted(p.name for p in here.iterdir() if p.is_dir()))
