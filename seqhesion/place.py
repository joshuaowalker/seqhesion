"""Placing ITS2-only sequences beside their closest tips (Josh, 2026-09-30).

ITS2-only inputs never shape the hierarchy; they are placed on it afterwards:
  * ITS2 is cut out of every tip (pyitsx extract --region ITS2; a tip is oriented full ITS, and its
    ITS2 matched the ITS2 cut from the raw input in 299 of 300 checked) and out of every query;
  * each query is compared with the tips' ITS2 (vsearch --usearch_global, global identity, --iddef 2);
    its closest tips are all those at the best identity (ties are common: ITS2 resolves less than
    full ITS);
  * it is placed in the finest group that holds every one of those closest tips. If they span
    groups even at the coarsest level, it is placed in their component only.
Reported per query: the closest tips' identity and how many there are, the group and level placed
at, and the margin to the best match outside that group.
"""
import collections
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from .fasta import read_fasta, write_fasta

PYITSX = os.path.join(os.path.dirname(sys.executable), 'pyitsx')


def extract_its2(seqs, cpus=8):
    """{id: ITS2} for those of `seqs` ({id: sequence}) in which pyitsx finds an ITS2."""
    with tempfile.TemporaryDirectory() as d:
        write_fasta(f'{d}/in.fasta', seqs)
        subprocess.run([PYITSX, 'extract', '-i', f'{d}/in.fasta', '-o', f'{d}/out.fasta', '--region', 'ITS2', '--cpus', str(cpus)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return {k.split('|')[0]: v.upper() for k, v in read_fasta(f'{d}/out.fasta').items() if v}


def closest(queries, refs, threads=8, exclude=None):
    """For each query ({id: seq}): (best identity, [ref ids at that identity], best identity among
    refs NOT in that list... computed by the caller). Refs ({id: seq}) are dereplicated first, so a
    tie between identical reference sequences costs one hit. `exclude`: {query: ref id} never to
    report (leave-one-out). Returns {query: [(identity, [ref ids]) ...] best first}."""
    by_seq = collections.defaultdict(list)
    for r, s in refs.items():
        by_seq[s].append(r)
    uniq = {f'r{k}': s for k, s in enumerate(by_seq)}
    members = {f'r{k}': ids for k, ids in enumerate(by_seq.values())}
    out = collections.defaultdict(list)
    with tempfile.TemporaryDirectory() as d:
        write_fasta(f'{d}/q.fasta', queries)
        write_fasta(f'{d}/db.fasta', uniq)
        subprocess.run(['vsearch', '--usearch_global', f'{d}/q.fasta', '--db', f'{d}/db.fasta', '--id', '0.70', '--iddef', '2',
                        '--maxaccepts', '64', '--maxrejects', '256', '--threads', str(threads),
                        '--userout', f'{d}/hits.tsv', '--userfields', 'query+target+id', '--quiet'],
                       stderr=subprocess.DEVNULL, check=True)
        for line in open(f'{d}/hits.tsv'):
            q, t, i = line.rstrip('\n').split('\t')
            ids = [r for r in members[t] if not exclude or r != exclude.get(q)]
            if ids:
                out[q].append((float(i) / 100, ids))
    for q in out:
        out[q].sort(key=lambda x: -x[0])
    return dict(out)


def place(hits, levels, group_at):
    """hits: closest() for one query; levels: finest first; group_at[level][tip] -> group id or None.
    Returns (identity, closest tips, level placed at or None, group or None, margin): margin is the
    best identity minus the best identity of any tip outside the placed group (None if none seen)."""
    if not hits:
        return None
    best = hits[0][0]
    top = [r for i, ids in hits if i == best for r in ids]
    level = group = None
    for L in levels:
        gs = {group_at[L].get(t) for t in top}
        if len(gs) == 1 and None not in gs:
            level, group = L, gs.pop()
            break
    outside = [i for i, ids in hits for r in ids if group is None and r not in top or group is not None and group_at[level].get(r) != group]
    margin = best - max(outside) if outside else None
    return best, sorted(top), level, group, margin
