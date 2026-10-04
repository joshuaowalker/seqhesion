"""Placing a sequence onto a component's built shard trees, without changing them.

Used for ITS2-only inputs (which never shape the hierarchy) and, later, for daily increments of new
full-ITS sequences. For a query q and a candidate tip c (found by an approximate search; identity
is used only to choose where to look):
  * shards: the SHARDS_PER_COVER shards of each cover that hold c, c most central first (the seed,
    then by neighbour rank);
  * per shard: q is aligned onto the shard's stored alignment (`mafft --add`, or `--addfragments`
    for a partial sequence such as ITS2 alone, with `--keeplength`), the shard's own column mask
    applied, and a FastTree tree re-inferred on the stored trimmed alignment plus q. q's join level
    with every member is read from that tree, prepared as the hierarchy reads its trees. The stored
    tree, and so every existing pair's evidence, is untouched;
  * pooled: q's join level with each member is the median over the shards holding both;
  * assigned: the finest level at which some group's average join level to q is within the level
    (the average linkage the hierarchy is built with), that group, and its ancestors above.
Diagnostic: whether the re-inference rearranged q's surroundings relative to the stored tree
(collapsed trees: q's attachment, an edge or a polytomy, does not exist in the stored tree).

Leave-one-out on 100 tips (lab, experiments/seqhesion/RESULTS.md, 2026-10-02): a tip's ITS2 alone
is assigned its own group 78-96% of the time across levels (vsearch identity placement: 53-97%);
its full ITS, 96-99%.
"""
import collections
import json
import os
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import coassoc, shardtrees, tools
from .fasta import read_fasta, write_fasta

SHARDS_PER_COVER = 3


def shard_mask(trees_dir, shard):
    """Columns of the stored alignment that trimAl -gappyout kept (the shard's trimmed alignment),
    cached beside it; checked once against the stored trimmed alignment."""
    td = Path(trees_dir)
    cache = td / f'{shard}.mask.json'
    if cache.exists():
        return json.load(open(cache))
    out = subprocess.run([tools.TRIMAL, '-in', str(td / f'{shard}.aln.fasta'), '-out', '/dev/null', '-gappyout', '-colnumbering'],
                         capture_output=True, text=True, check=True).stdout
    cols = [int(x) for x in out.split('\t', 1)[1].replace(' ', '').strip().split(',')]
    aln, trim = read_fasta(td / f'{shard}.aln.fasta'), read_fasta(td / f'{shard}.trim.fasta')
    if any(''.join(aln[k][c] for c in cols) != trim[k] for k in trim):
        raise RuntimeError(f'{shard}: the recomputed trimAl mask does not reproduce the stored trimmed alignment')
    tmp = Path(f'{cache}.{os.getpid()}.tmp')          # workers may compute the same shard's mask at once
    json.dump(cols, open(tmp, 'w'))
    tmp.replace(cache)
    return cols


def central_shards(manifest, tip, per_cover=SHARDS_PER_COVER):
    """The shards of each cover holding `tip`, where it is most central (its rank among the members:
    the seed first, then its neighbours by identity). A cover-1 shard that is a twin of a cover-0
    shard (cover.twins: the same sample) is skipped."""
    from .cover import twins
    skip = twins(manifest)
    out = []
    for c in (0, 1):
        held = sorted((s['members'].index(tip), name) for name, s in manifest['shards'].items()
                      if s['cover'] == c and tip in s['members'] and name not in skip)
        out += [name for _, name in held[:per_cover]]
    return out


def leafset(n):
    return frozenset(l.name for l in n.iter_leaves())


def collapse(t, min_len):
    for n in list(t.traverse('postorder')):
        if not n.is_leaf() and not n.is_root() and n.dist < min_len:
            n.delete(prevent_nondicotomic=False, preserve_branch_length=True)
    return t


def _sides(n, everything):
    sides = {leafset(c) for c in n.children}
    if not n.is_root():
        sides.add(everything - leafset(n))
    return sides


def rearranged(stored_path, reinferred_path, q, min_len):
    """True if q's attachment in the re-inference (collapsed) has no counterpart in the stored tree."""
    from ete3 import Tree
    R = collapse(Tree(str(reinferred_path)), min_len)
    T = collapse(Tree(str(stored_path)), min_len)
    everything = leafset(T)
    P = (R & q).up
    want = {s - {q} for s in _sides(P, leafset(R))} - {frozenset()}
    want = {s & everything for s in want} - {frozenset()}
    if any(_sides(n, everything) == want for n in T.traverse()):
        return False
    if len(want) == 2:
        U, V = want
        return not any(not n.is_root() and leafset(n) in (U, V) for n in T.traverse())
    return True


def join_values(tree_path, trim_path, members, scaffold, q):
    """{member: q's join level with it} in the prepared tree ('identical' join)."""
    _, line, _, _ = shardtrees.prepare(('q', tree_path, trim_path, members, scaffold))
    names = sorted(members)
    index = {t: i for i, t in enumerate(names)}
    coassoc._init(index, 'identical')
    coassoc._G['floor'] = 1.0 / tools.alignment_width(trim_path)
    keys, vals = coassoc.joins(line)
    n, qi = len(index), index[q]
    out = {}
    for k, v in zip(keys.tolist(), vals.tolist()):
        i, j = divmod(k, n)
        if i == qi:
            out[names[j]] = v
        elif j == qi:
            out[names[i]] = v
    return out


def in_shard(job):
    """q's join levels in one shard, and whether the re-inference rearranged its surroundings."""
    trees_dir, shard, q, seq, fragment, members, scaffold = job
    td = Path(trees_dir)
    cols = shard_mask(td, shard)
    trim = read_fasta(td / f'{shard}.trim.fasta')
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        write_fasta(d / 'q.fasta', {q: seq})
        with open(d / 'add.aln', 'w') as fh:
            subprocess.run([tools.MAFFT, '--addfragments' if fragment else '--add', str(d / 'q.fasta'), '--keeplength',
                            '--quiet', '--thread', '1', str(td / f'{shard}.aln.fasta')],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True, env=tools._mafft_env())
        row = read_fasta(d / 'add.aln')[q]
        write_fasta(d / 'all.trim', {**trim, q: ''.join(row[c] for c in cols)})
        tools.run_fasttree(d / 'all.trim', d / 'q.nwk')
        min_len = 0.5 / len(next(iter(trim.values())))
        moved = rearranged(td / f'{shard}.nwk', d / 'q.nwk', q, min_len)
        joins = join_values(d / 'q.nwk', td / f'{shard}.trim.fasta', list(members) + [q], scaffold, q)
    return q, shard, joins, moved


def assign(pooled, group_at, levels):
    """[(group or None) per level]: the finest level at which a group's average join level to q is
    within the level, that group, and its ancestors above."""
    for k, L in enumerate(levels):
        by = collections.defaultdict(list)
        for y, v in pooled.items():
            g = group_at[L].get(y)
            if g is not None:
                by[g].append(v)
        if not by:
            continue
        g, avg = min(((g, float(np.mean(v))) for g, v in by.items()), key=lambda t: (t[1], t[0]))
        if avg <= L:
            tips = {t for t, h in group_at[L].items() if h == g}
            out = [None] * k
            for L2 in levels[k:]:
                anc = {group_at[L2].get(t) for t in tips}
                out.append(anc.pop() if len(anc) == 1 else None)
            return out, avg
    return [None] * len(levels), None


def place(queries, covers_dir, group_at, levels, procs=8, log=print):
    """queries: [(id, sequence, fragment, candidate tip)] for one component. Returns
    {id: {'groups': [per level], 'join': average join level to the finest group, 'shards': n,
    'rearranged': n}}."""
    covers_dir = Path(covers_dir)
    man = json.load(open(covers_dir / 'manifest.json'))
    jobs = []
    for q, seq, fragment, cand in queries:
        for shard in central_shards(man, cand):
            jobs.append((str(covers_dir / 'trees'), shard, q, seq, fragment, man['shards'][shard]['members'], man['scaffold']))
    pooled = collections.defaultdict(lambda: collections.defaultdict(list))
    shards, moved = collections.Counter(), collections.Counter()
    with Pool(procs) as pool:
        for i, (q, shard, joins, mv) in enumerate(pool.imap_unordered(in_shard, jobs, chunksize=1), 1):
            shards[q] += 1
            moved[q] += mv
            for y, v in joins.items():
                pooled[q][y].append(v)
            if i % 500 == 0:
                log(f'  {i}/{len(jobs)} shard placements')
    out = {}
    for q, _, _, _ in queries:
        med = {y: float(np.median(v)) for y, v in pooled[q].items()}
        groups, avg = assign(med, group_at, levels)
        out[q] = {'groups': groups, 'join': avg, 'shards': shards[q], 'rearranged': moved[q]}
    return out
