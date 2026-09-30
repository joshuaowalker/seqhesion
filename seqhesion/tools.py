"""One shard tree: MAFFT L-INS-i -> trimAl -gappyout -> FastTree -gtr -gamma, cached.

The recipe and its cache key are carried over unchanged from the ubertree lab, so trees built
there are reused here and the reverse; `build_key` must fingerprint exactly as it did there
(experiments/stitch/RESULTS.md in the lab has the verification of the MAFFT build).
"""
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

from . import cache
from .fasta import read_fasta, write_fasta

# Prefer the Apple-silicon build of MAFFT when it is installed: upstream 7.526 plus bit-exact
# rewrites, verified byte-identical to stock on ~4,500 windows for this invocation (L-INS-i)
# only. Stock and this build share one tree cache, on purpose: the key uses the upstream version
# (mafft_key_version). The build began reporting `v7.526-opt4` on 2026-09-28; the 12-window canary
# (lab: experiments/stitch/verify_mafft.py --largest) was byte-identical on 2026-09-30, so the
# suffix names a build, not a different result. If a build ever diverges, key on the suffix and
# rebuild on purpose. ~/mm/local/bin is deliberately not on PATH.
_MAFFT_FAST = os.path.expanduser('~/mm/local/bin/mafft')
MAFFT = _MAFFT_FAST if os.access(_MAFFT_FAST, os.X_OK) else (shutil.which('mafft') or '/opt/homebrew/bin/mafft')
TRIMAL = shutil.which('trimal') or '/opt/homebrew/bin/trimal'
_CONDA_BIN = os.path.expanduser('~/mm/analysis/phylogeny/conda-env/bin')
FASTTREE = shutil.which('FastTree') or f'{_CONDA_BIN}/FastTree'

# Recorded on every tree as provenance, not part of the key (an OS upgrade was verified not to
# change CPU MAFFT output; keying on it would invalidate every tree on every upgrade).
_OS_BUILD = platform.platform()

# Bump when the commands run by build_tree change in a way that changes results.
RECIPE = 'mafft --localpair --maxiterate 1000 | trimal | FastTree -gtr -gamma -nt; v1'


def _mafft_env():
    """MAFFT_BINARIES would make mafft run other binaries than the verified build."""
    env = dict(os.environ)
    env.pop('MAFFT_BINARIES', None)
    return env


def run_mafft(in_fasta, out_fasta, threads=1):
    """L-INS-i as MAFFT defines it. Refinement converges after 4-12 cycles on our windows, so
    any --maxiterate >= 16 gives identical output. Returns stderr, which under --quiet is empty
    for a normal run, so anything there is signal."""
    with open(out_fasta, 'wb') as out:
        p = subprocess.run([MAFFT, '--localpair', '--maxiterate', '1000', '--quiet',
                            '--thread', str(threads), str(in_fasta)],
                           stdout=out, stderr=subprocess.PIPE, check=True, env=_mafft_env())
    return p.stderr.decode('utf-8', 'replace').strip()


def run_trimal(in_fasta, out_fasta):
    subprocess.run([TRIMAL, '-in', str(in_fasta), '-out', str(out_fasta), '-gappyout'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    write_fasta(out_fasta, read_fasta(out_fasta))      # trimAl appends " NNN bp" to headers; keep bare ids


def run_fasttree(in_fasta, out_tree):
    with open(out_tree, 'wb') as out, open(os.devnull, 'w') as err:
        subprocess.run([FASTTREE, '-gtr', '-gamma', '-nt', '-quiet', str(in_fasta)], stdout=out, stderr=err, check=True)


def alignment_width(path):
    return len(next(iter(read_fasta(path).values())))


def build_tree(seqs, workdir, name, order=None):
    """Infer a tree for `seqs` ({id: sequence}) in `order`. Writes <name>.fasta, .aln.fasta,
    .trim.fasta and .nwk under `workdir`. Returns (tree_path, info)."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    raw, aln, trimmed, tree = (workdir / f'{name}.{e}' for e in ('fasta', 'aln.fasta', 'trim.fasta', 'nwk'))
    write_fasta(raw, seqs, order)
    info = {'n': len(seqs)}
    warned = run_mafft(raw, aln)
    if warned:
        info['aligner_stderr'] = warned[:2000]
    info['cols_aligned'] = alignment_width(aln)
    run_trimal(aln, trimmed)
    info['cols_trimmed'] = alignment_width(trimmed)
    tmp = workdir / f'{name}.nwk.tmp'           # never leave a truncated tree under the final name
    run_fasttree(trimmed, tmp)
    os.replace(tmp, tree)
    return tree, info


def mafft_key_version():
    """The MAFFT version as the cache key sees it: upstream's, without a verified build's suffix."""
    return re.sub(r'-opt\d+$', '', cache.tool_version(MAFFT))


def build_key(seqs, order=None):
    """Fingerprint of everything that determines build_tree's output. The fields `trim`,
    `aligned` and `tree_method` are fixed here but kept in the key, so that it matches the keys
    of trees the lab already built."""
    tools = {'mafft': mafft_key_version(), 'trimal': cache.tool_version(TRIMAL),
             'tree': cache.tool_version(FASTTREE)}
    return cache.fingerprint(sequences=cache.sequences_hash(seqs, order), trim=True, aligned=False,
                             tree_method='fasttree', recipe=RECIPE, tools=tools)


def cached_build_tree(seqs, workdir, name, order=None):
    """build_tree, reusing <workdir>/<name>.nwk only if its key says it was made from exactly
    these sequences, in this order, with this recipe and these tools (cache.Stale otherwise).
    Returns (tree_path, reused)."""
    workdir = Path(workdir)
    tree = workdir / f'{name}.nwk'
    key = build_key(seqs, order)
    ids = list(order or seqs)
    if cache.valid(tree, key, tips=ids):
        return tree, True
    _, info = build_tree(seqs, workdir, name, order=order)
    cache.record(tree, key, {'n': len(ids), 'recipe': RECIPE, 'aligner': 'linsi', 'aligner_binary': MAFFT,
                             'aligner_version': cache.tool_version(MAFFT),
                             **({'aligner_stderr': info['aligner_stderr']} if 'aligner_stderr' in info else {}),
                             'os_build': _OS_BUILD})
    return tree, False
