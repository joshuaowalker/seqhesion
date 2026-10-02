"""Placing ITS2-only inputs on the hierarchy (Josh, 2026-09-30 / 2026-10-02).

ITS2-only inputs never shape the hierarchy; they are placed on it afterwards, through the shard
trees (seqhesion.insert), never by an identity threshold:
  * ITS2 is cut out of every query and every tip (pyitsx extract --region ITS2);
  * candidate search (sampling only): the query's ITS2 against all tips' ITS2 (vsearch); the
    candidate is the closest tip (ties: the first id);
  * the query's ITS2 is added as a fragment to the shards where the candidate is central, and
    assigned to the finest group whose average join level to it is within the level (and that
    group's ancestors).
Statuses: placed (in a group at some level) / ungrouped (in the component, but no group within any
level) / not_built (the candidate's component is not in this release) / spans_components (equally
close candidates in more than one component) / no_match (no tip ITS2 within the search's 70%
identity floor) / no_its2 (pyitsx found no ITS2 in the query).
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


PLACE_COLUMNS = ['input_id', 'status', 'component', 'group_id', 'level', 'join', 'shards', 'rearranged',
                 'candidate', 'candidate_identity']


def cached_its2(seqs, path, cpus=8):
    """extract_its2 with its result kept at `path` (a FASTA) beside the intake; reused only if it
    covers exactly the requested ids (ids are content hashes, so same id = same sequence)."""
    path = Path(path)
    if path.exists():
        have = read_fasta(path)
        done = set(read_fasta(f'{path}.ids')) if Path(f'{path}.ids').exists() else set()
        if done == set(seqs):
            return have
    got = extract_its2(seqs, cpus)
    write_fasta(f'{path}.tmp', got, sorted(got))
    write_fasta(f'{path}.ids.tmp', {i: 'N' for i in seqs}, sorted(seqs))
    os.replace(f'{path}.tmp', path)
    os.replace(f'{path}.ids.tmp', f'{path}.ids')
    return got


def place_inputs(intake_dir, classes, built_group_at, levels, comp_of, covers_of, threads=8, log=print):
    """Rows of placements.tsv for every ITS2-only input. covers_of: {component: covers directory}
    for the components this release built. The candidate search compares with ALL tips of the
    intake, so a query whose relatives lie in an unbuilt component is reported as not_built rather
    than placed beside a distant built tip."""
    from . import insert
    intake_dir = Path(intake_dir)
    tips = read_fasta(intake_dir / 'full.derep.fasta')
    raw = read_fasta(intake_dir / 'readable.fasta')
    queries = {r['seq_id']: raw[r['seq_id']] for r in classes if r['class'] == 'its2'}
    tip_its2 = cached_its2(tips, intake_dir / 'tips.its2.fasta', threads)
    q_its2 = cached_its2(queries, intake_dir / 'its2_inputs.its2.fasta', threads)
    hits = closest(q_its2, tip_its2, threads)
    rows = {}
    todo = collections.defaultdict(list)
    for q in sorted(queries):
        row = dict.fromkeys(PLACE_COLUMNS)
        row['input_id'] = q
        rows[q] = row
        if q not in q_its2:
            row['status'] = 'no_its2'
            continue
        if q not in hits:
            row['status'] = 'no_match'
            continue
        best = hits[q][0][0]
        top = sorted({r for i, ids in hits[q] if i == best for r in ids})
        comps = {comp_of.get(t) for t in top}
        row.update(candidate=top[0], candidate_identity=round(best, 4))
        if len(comps) > 1:
            row['status'] = 'spans_components'
            continue
        comp = comps.pop()
        row['component'] = comp
        if comp not in covers_of:
            row['status'] = 'not_built'
            continue
        todo[comp].append((q, q_its2[q], True, top[0]))
    for comp, qs in todo.items():
        log(f'placing {len(qs)} ITS2-only inputs in component {comp}')
        res = insert.place(qs, covers_of[comp], built_group_at, levels, threads, log)
        for q, r in res.items():
            k = next((i for i, g in enumerate(r['groups']) if g is not None), None)
            rows[q].update(status='placed' if k is not None else 'ungrouped',
                           group_id=r['groups'][k] if k is not None else None, level=levels[k] if k is not None else None,
                           join=None if r['join'] is None else round(r['join'], 5), shards=r['shards'], rearranged=r['rearranged'])
    status = collections.Counter(r['status'] for r in rows.values())
    log(f'ITS2-only placements: {dict(status)}')
    return [rows[q] for q in sorted(rows)], dict(status)
