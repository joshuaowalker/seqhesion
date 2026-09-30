"""Reusing results safely.

A product on disk (a tree, an alignment) may be reused only if it was made from the
same inputs, by the same steps, with the same tools. Every product gets a sidecar
`<product>.key.json` holding a fingerprint of all three; a product whose sidecar
does not match is *stale* and reusing it raises instead of silently handing back a
tree of other sequences (file names say nothing: shard `c0.0042` means different
members after any change to the cover builder).

Products are written to a temporary name and renamed at the end, so an interrupted
job never leaves a truncated file that looks finished.
"""
import functools
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


class Stale(Exception):
    """A product exists but was not made from what is being asked for now."""


def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(**parts):
    """sha256 over a canonical JSON rendering; Path values are hashed by content."""
    def norm(v):
        if isinstance(v, Path):
            return {'file': file_hash(v)}
        if isinstance(v, dict):
            return {str(k): norm(x) for k, x in sorted(v.items())}
        if isinstance(v, (list, tuple)):
            return [norm(x) for x in v]
        return v
    return hashlib.sha256(json.dumps(norm(parts), sort_keys=True).encode()).hexdigest()


def sequences_hash(seqs, order=None):
    """Identity of an input set: ids and sequences in the order they are written
    (MAFFT's result depends on input order)."""
    h = hashlib.sha256()
    for i in (order or seqs):
        h.update(f'>{i}\n{seqs[i].upper()}\n'.encode())
    return h.hexdigest()


@functools.lru_cache(maxsize=None)
def tool_version(path):
    for flag in ('--version', '-version', '-expert'):
        try:
            p = subprocess.run([path, flag], capture_output=True, text=True, timeout=20, stdin=subprocess.DEVNULL)
        except (OSError, subprocess.TimeoutExpired):
            continue
        m = re.search(r'(?:version|\bv|FastTree)\s*([0-9][0-9A-Za-z.\-]*( Double precision)?)', p.stdout + p.stderr, re.I)
        if m:
            return m.group(1)
    return 'unknown'


def sidecar(product):
    return Path(f'{product}.key.json')


def newick_tips(path):
    return set(re.findall(r'[(,]([^(),:;]+):', Path(path).read_text()))


def valid(product, key, tips=None):
    """True: reuse it. False: not there, build it. Raises Stale when something is
    there that was made from different inputs, or has no record of what made it."""
    product = Path(product)
    if not product.exists():
        return False
    side = sidecar(product)
    if not side.exists():
        raise Stale(f'{product} exists without a key (made before keys were recorded); adopt or remove it')
    have = json.load(open(side))['key']
    if have != key:
        raise Stale(f'{product} was made from different inputs or settings (key {have[:12]}, wanted {key[:12]})')
    if tips is not None and newick_tips(product) != set(tips):
        raise Stale(f'{product} has the right key but its tips differ from the request')
    return True


def commit(tmp, product, key, meta=None):
    """Move a finished temporary file into place and record what made it."""
    os.replace(tmp, product)
    record(product, key, meta)


def record(product, key, meta=None):
    tmp = f'{sidecar(product)}.tmp'
    json.dump({'key': key, **(meta or {})}, open(tmp, 'w'), indent=1, sort_keys=True)
    os.replace(tmp, sidecar(product))
