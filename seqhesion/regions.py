"""Splitting the corpus into regions: independent pieces, each built on its own.

A region is a connected component of the centroid graph: 90% centroids of the dereplicated
full-ITS set (`vsearch --cluster_fast`, used only to locate regions, never structurally), linked
when two centroids match at >= min_id. Nothing is chosen or dropped by name.

A region directory holds:
  comp.fasta             its sequences, most copies first
  knn.tsv                vsearch --usearch_global --self --id 0.80: the neighbourhoods shards are made of
  comp_cent_0.85.fasta   85% centroids with sizes: scaffold candidates (rooting context)
  region.json            how it was made
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


def components(intake_dir, min_id=80.0):
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


def regions(intake_dir, out_dir, min_id=80.0, min_size=2, threads=6, log=print):
    """Every region of 2+ sequences under out_dir/r0000, r0001, ... (largest first)."""
    centroid_graph(intake_dir, threads)
    seqs = read_fasta(Path(intake_dir) / 'full.derep.fasta')
    sizes = copies(intake_dir)
    comps = components(intake_dir, min_id)
    kept = [c for c in comps if len(c) >= min_size]
    log(f'{len(comps)} regions; {len(kept)} of {min_size}+ sequences hold {sum(map(len, kept))} of {len(seqs)}')
    for k, ids in enumerate(kept):
        make_region({i: seqs[i] for i in ids}, sizes, Path(out_dir) / f'r{k:04d}', threads,
                    how={'intake': str(intake_dir), 'min_id': min_id, 'component': k})
    return kept
