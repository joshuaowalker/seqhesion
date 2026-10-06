"""Splitting the corpus into regions: independent pieces, each built on its own.

A region is a connected component of the centroid graph: 90% centroids of the dereplicated
full-ITS set (`vsearch --cluster_fast`, used only to locate regions, never structurally), linked
when two centroids match at >= min_id (default 86: at 80, most of the Agaricales falls into one
component of ~60K sequences). Nothing is chosen or dropped by name.

Region directories are named by their anchor (the tip with most copies), so a rebuild of a grown
corpus finds the directories, and the cached trees, of components it has built before. A component
of SHARD or more tips is LARGE; one of 2 to SHARD - 1 tips is SMALL (below); a lone tip gets no
region.

A region directory holds:
  comp.fasta             its sequences, most copies first
  knn.tsv                vsearch --usearch_global --self --id 0.80: the neighbourhoods shards are made of
  comp_cent_0.85.fasta   85% centroids with sizes: scaffold candidates (rooting context)
  region.json            how it was made

A SMALL region (a component smaller than one shard; decided 2026-10-03) also holds
  context.fasta          the outside sequences its shards draw on
and its knn.tsv lists each tip's neighbours in the WHOLE corpus (vsearch floor SMALL_KNN_ID), so its
shards are seeded from its own tips but filled with their nearest sequences anywhere. Those outside
sequences are rooting context only: pruned before co-association, never grouped. Built alone, such a
component's shards hold all of it and its only rooting candidates are its own tips (lab
experiment, 2026-10-03: 0% of its tips join another component's tips at levels
<= 0.075, so pruning loses nothing). Its scaffold is the neighbourhood's 85% centroids.
"""
import collections
import json
import os
import re
import subprocess
from pathlib import Path

from . import cache
from .fasta import read_fasta, write_fasta

KNN = ['--id', '0.80', '--maxaccepts', '400', '--maxrejects', '64']
SMALL_KNN_ID = 0.60
SHARD = 150                     # cover.plan's shard size: a component smaller than this is small
QUOTA = 30


def _strip(h):
    return re.sub(r';size=\d+;?$', '', h)


def centroid_graph(intake_dir, threads=6):
    """clust_0.90.uc, cent_0.90.fasta and cent90_graph.tsv (centroids against centroids,
    --id 0.70, up to 40 hits each) beside the intake's full.derep.fasta."""
    d = Path(intake_dir)
    derep = d / 'full.derep.fasta'
    cent, uc, graph = d / 'cent_0.90.fasta', d / 'clust_0.90.uc', d / 'cent90_graph.tsv'
    if not uc.exists():
        subprocess.run(['vsearch', '--cluster_fast', str(derep), '--id', '0.90', '--centroids', f'{cent}.tmp', '--uc', f'{uc}.tmp',
                        '--threads', str(threads), '--quiet'], check=True)
        os.replace(f'{cent}.tmp', cent)
        os.replace(f'{uc}.tmp', uc)
    if not graph.exists():
        subprocess.run(['vsearch', '--usearch_global', str(cent), '--db', str(cent), '--self', '--id', '0.70', '--maxaccepts', '40',
                        '--maxrejects', '64', '--threads', str(threads), '--userout', f'{graph}.tmp',
                        '--userfields', 'query+target+id', '--quiet'], check=True)
        os.replace(f'{graph}.tmp', graph)


MIN_ID = 86.0


def components(intake_dir, min_id=MIN_ID):
    """Every region: lists of sequence ids, largest first."""
    d = Path(intake_dir)
    members = collections.defaultdict(list)                     # centroid -> sequences
    for line in open(d / 'clust_0.90.uc'):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'S':
            members[_strip(f[8])].append(_strip(f[8]))
        elif f[0] == 'H':
            members[_strip(f[9])].append(_strip(f[8]))
    parent = {c: c for c in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for line in open(d / 'cent90_graph.tsv'):
        q, t, ident = line.rstrip('\n').split('\t')
        if float(ident) >= min_id:
            q, t = _strip(q), _strip(t)
            if q in parent and t in parent:
                parent[find(q)] = find(t)
    comps = collections.defaultdict(list)
    for c in members:
        comps[find(c)].extend(members[c])
    return sorted((sorted(v) for v in comps.values()), key=lambda v: (-len(v), v[0]))


def copies(intake_dir):
    """{unique id: input records carrying it} from the intake's derep_members.tsv."""
    out = {}
    for line in open(Path(intake_dir) / 'derep_members.tsv'):
        if not line.startswith('unique\t'):
            u, size, _ = line.rstrip('\n').split('\t', 2)
            out[u] = int(size)
    return out


def make_region(seqs, sizes, out, threads=6, how=None):
    """Write a region directory for `seqs` ({id: sequence}), `sizes` ({id: copies})."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    order = sorted(seqs, key=lambda i: (-sizes[i], i))
    write_fasta(out / 'comp.fasta', seqs, order)
    subprocess.run(['vsearch', '--usearch_global', str(out / 'comp.fasta'), '--db', str(out / 'comp.fasta'), '--self', *KNN,
                    '--threads', str(threads), '--userout', str(out / 'knn.tsv'), '--userfields', 'query+target+id', '--quiet'],
                   stderr=subprocess.DEVNULL, check=True)
    sized = out / 'comp.sized.fasta'
    write_fasta(sized, {f'{i};size={sizes[i]}': seqs[i] for i in order})
    subprocess.run(['vsearch', '--cluster_fast', str(sized), '--id', '0.85', '--centroids', str(out / 'comp_cent_0.85.fasta'),
                    '--sizein', '--sizeout', '--threads', str(threads), '--quiet'], stderr=subprocess.DEVNULL, check=True)
    sized.unlink()
    json.dump({**(how or {}), 'sequences': len(seqs), 'sequences_hash': cache.sequences_hash(seqs, order)[:16],
               'knn': 'vsearch --usearch_global ' + ' '.join(KNN)}, open(out / 'region.json', 'w'), indent=1)


def regions(intake_dir, out_dir, min_id=MIN_ID, threads=6, log=print):
    """A region dir under out_dir for every component of 2+ tips, named by its anchor; existing
    ones are kept. Writes out_dir/large.txt and small.txt (region dirs, largest component first)
    and returns (large, small)."""
    intake_dir, out_dir = Path(intake_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    centroid_graph(intake_dir, threads)
    seqs = read_fasta(intake_dir / 'full.derep.fasta')
    sizes = copies(intake_dir)
    comps = components(intake_dir, min_id)                  # largest first
    big = [c for c in comps if len(c) >= SHARD]
    log(f'{len(comps)} components: {len(big)} large ({sum(map(len, big))} tips), '
        f'{sum(1 for c in comps if 2 <= len(c) < SHARD)} small, {sum(1 for c in comps if len(c) == 1)} lone tips')
    large = []
    for ids in big:
        rd = out_dir / min(ids, key=lambda t: (-sizes[t], t))
        if not (rd / 'region.json').exists():
            make_region({t: seqs[t] for t in ids}, sizes, rd, threads, how={'intake': str(intake_dir), 'min_id': min_id})
            log(f'  made {rd.name}: {len(ids)} tips')
        large.append(rd)
    small = small_regions(intake_dir, comps, out_dir, threads, log=log)
    for name, dirs in (('large.txt', large), ('small.txt', small)):
        (out_dir / name).write_text(''.join(f'{d}\n' for d in dirs))
    return large, small


def small_regions(intake_dir, comps, out_dir, threads=6, log=print):
    """Region dirs for small components (2 to SHARD - 1 tips): one corpus-wide neighbour search for
    all their tips, then per component comp.fasta, context.fasta, knn.tsv, comp_cent_0.85.fasta and
    region.json (kind 'small'). Directories are named by anchor, the tip with most copies. Returns them."""
    d, out_dir = Path(intake_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(d / 'full.derep.fasta')
    sizes = copies(intake_dir)
    comps = [sorted(c) for c in comps if 2 <= len(c) < SHARD]
    dirs = [out_dir / min(c, key=lambda t: (-sizes[t], t)) for c in comps]
    todo = [(c, rd) for c, rd in zip(comps, dirs) if not (rd / 'region.json').exists()]
    if not todo:
        return dirs
    q, hits = out_dir / 'small_queries.fasta', out_dir / 'small_knn_corpus.tsv'
    write_fasta(q, seqs, sorted(t for c, _ in todo for t in c))
    subprocess.run(['vsearch', '--usearch_global', str(q), '--db', str(d / 'full.derep.fasta'), '--self', '--id', str(SMALL_KNN_ID),
                    '--maxaccepts', '400', '--maxrejects', '64', '--threads', str(threads), '--userout', f'{hits}.tmp',
                    '--userfields', 'query+target+id', '--quiet'], stderr=subprocess.DEVNULL, check=True)
    os.replace(f'{hits}.tmp', hits)
    rows = collections.defaultdict(list)
    for line in open(hits):
        rows[line.split('\t', 1)[0]].append(line)
    for c, rd in todo:
        rd.mkdir(parents=True, exist_ok=True)
        C = set(c)
        with open(rd / 'knn.tsv', 'w') as f:
            for t in c:
                f.writelines(rows.get(t, []))
        ranked = {t: [x for _, x in sorted(((float(l.rstrip().split('\t')[2]), l.split('\t')[1]) for l in rows.get(t, [])),
                                          reverse=True)] for t in c}
        nbhd = C | {x for t in c for x in ranked[t][:SHARD - 1 + QUOTA]}
        outside = sorted(nbhd - C)
        write_fasta(rd / 'comp.fasta', {t: seqs[t] for t in c}, sorted(c, key=lambda i: (-sizes[i], i)))
        write_fasta(rd / 'context.fasta', {t: seqs[t] for t in outside}, outside)
        sized = rd / 'nbhd.sized.fasta'
        order = sorted(nbhd, key=lambda i: (-sizes.get(i, 1), i))
        write_fasta(sized, {f'{i};size={sizes.get(i, 1)}': seqs[i] for i in order})
        subprocess.run(['vsearch', '--cluster_fast', str(sized), '--id', '0.85', '--centroids', str(rd / 'comp_cent_0.85.fasta'),
                        '--sizein', '--sizeout', '--threads', str(threads), '--quiet'], stderr=subprocess.DEVNULL, check=True)
        sized.unlink()
        json.dump({'kind': 'small', 'intake': str(intake_dir), 'sequences': len(c), 'context': len(outside),
                   'sequences_hash': cache.sequences_hash({t: seqs[t] for t in c}, sorted(c, key=lambda i: (-sizes[i], i)))[:16],
                   'knn': f'vsearch --usearch_global against the whole corpus --id {SMALL_KNN_ID} --maxaccepts 400 --maxrejects 64'},
                  open(rd / 'region.json', 'w'), indent=1)
    q.unlink()
    log(f'small regions: {len(todo)} made, {len(dirs) - len(todo)} already there')
    return dirs
