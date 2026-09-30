"""Fingerprints of what a hierarchy was built from, so a stale one is refused rather than believed.

A key is a dict of component digests, not one hash, so a mismatch can say which input moved.
Trees are identified by the sidecar `cache` wrote for each (content hash as a fallback).
"""
import json
from pathlib import Path

from . import cache

VERSION = 'derived-v1 (2026-09-28)'


def tree_key(path):
    side = cache.sidecar(path)
    if side.exists():
        return json.load(open(side))['key']
    return cache.file_hash(path)


def trees_key(paths):
    paths = sorted(Path(p) for p in paths)
    return cache.fingerprint(n=len(paths), trees={p.name: tree_key(p) for p in paths})


def hierarchy_key(covers_dir, method, settings):
    """The inputs of a hierarchy: the cover manifest, the shard trees of both covers, the prepared
    trees actually read (shardtrees.export), and the method and its settings."""
    covers_dir = Path(covers_dir)
    man = json.load(open(covers_dir / 'manifest.json'))
    trees = [covers_dir / 'trees' / f'{n}.nwk' for n in man['shards']]
    return {
        'version': VERSION,
        'manifest': cache.file_hash(covers_dir / 'manifest.json'),
        'trees': trees_key(trees),
        'exported': {f'cover{c}': cache.file_hash(covers_dir / 'coassoc' / f'cover{c}.nwk')[:16] for c in (0, 1)},
        'method': cache.fingerprint(method=method, **settings),
    }
