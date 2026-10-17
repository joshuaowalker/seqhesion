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
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

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


def components_digest(comp_of):
    """A fingerprint of the tip -> component map, so a placement plan is used only with the release
    whose components it was planned against."""
    import hashlib
    h = hashlib.sha256()
    for t in sorted(comp_of):
        h.update(f'{t}\t{comp_of[t]}\n'.encode())
    return h.hexdigest()


def plan_inputs(intake_dir, classes, comp_of, covers_of, threads=8, log=print):
    """Placement, step 1 (where the intake is): ITS2 cut from every query and tip, the candidate
    search, every status that needs no tree (no_its2, no_match, spans_components, not_built), and one
    job per (query, shard) for the rest. covers_of: {component: covers directory} for the components
    this release builds. The candidate search compares with ALL tips of the intake, so a query whose
    relatives lie in an unbuilt component is reported as not_built rather than placed beside a
    distant built tip."""
    from .insert import central_shards
    intake_dir = Path(intake_dir)
    tips = read_fasta(intake_dir / 'full.derep.fasta')
    raw = read_fasta(intake_dir / 'readable.fasta')
    queries = {r['seq_id']: raw[r['seq_id']] for r in classes if r['class'] == 'its2'}
    tip_its2 = cached_its2(tips, intake_dir / 'tips.its2.fasta', threads)
    q_its2 = cached_its2(queries, intake_dir / 'its2_inputs.its2.fasta', threads)
    hits = closest(q_its2, tip_its2, threads)
    rows, todo = {}, collections.defaultdict(list)
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
        todo[comp].append((q, top[0]))
    jobs, its2 = [], {}
    for comp, qs in todo.items():
        man = json.load(open(Path(covers_of[comp]) / 'manifest.json'))
        for q, cand in qs:
            its2[q] = q_its2[q]
            jobs += [[q, comp, shard] for shard in central_shards(man, cand)]
    jobs.sort(key=lambda j: (j[1], j[2], j[0]))       # by shard: a chunk touches few shards
    log(f'placement plan: {len(its2)} ITS2-only inputs to place by {len(jobs)} shard re-inferences in {len(todo)} '
        f'components; decided without trees: {dict(collections.Counter(r["status"] for r in rows.values() if r["status"]))}')
    return {'method': PLACE_METHOD, 'intake': str(intake_dir), 'components_sha256': components_digest(comp_of),
            'built': sorted(covers_of), 'covers': {c: str(covers_of[c]) for c in sorted(todo)},
            'rows': [rows[q] for q in sorted(rows)], 'its2': its2, 'jobs': jobs}


PLACE_METHOD = 'seqhesion.place v1: ITS2 fragment added to the 3 most central shards per cover of the closest tip'


def chunk_jobs(jobs, chunk, of):
    """Chunk `chunk` of `of`: a contiguous run of the (shard-sorted) jobs; the chunks partition them."""
    n = len(jobs)
    return jobs[chunk * n // of:(chunk + 1) * n // of]


def covers_dir(plan, comp, regions_root=None):
    """The covers directory of a component: as planned, or the same region and covers names under
    regions_root (a remote copy of the region dirs)."""
    p = Path(plan['covers'][comp])
    return Path(regions_root) / p.parent.name / p.name if regions_root else p


def run_jobs(plan, jobs, procs=8, regions_root=None, mask_cache=True, log=print):
    """Placement, step 2: every job's join levels; yields {'q', 'shard', 'joins', 'moved'} in job order."""
    from multiprocessing import Pool
    from . import insert
    mans, work = {}, []
    for q, comp, shard in jobs:
        cd = covers_dir(plan, comp, regions_root)
        if comp not in mans:
            mans[comp] = json.load(open(cd / 'manifest.json'))
        m = mans[comp]
        work.append((str(cd / 'trees'), shard, q, plan['its2'][q], True, m['shards'][shard]['members'], m['scaffold'], mask_cache))
    with Pool(procs) as pool:
        for i, (q, shard, joins, moved) in enumerate(pool.imap(insert.in_shard, work, chunksize=1), 1):
            if i % 500 == 0:
                log(f'  {i}/{len(work)} shard placements')
            yield {'q': q, 'shard': shard, 'joins': joins, 'moved': bool(moved)}


def place_shards(place_dir, regions_root, chunk, of, procs=8, log=print):
    """One chunk of a placement plan (chunks are independent: an AWS array): PLACE_DIR/joins/
    <chunk>-of-<of>.jsonl, one result per job. Needs only plan.json and, under REGIONS_ROOT, each
    planned shard's .aln.fasta / .trim.fasta / .nwk and its region's covers/manifest.json; writes
    nothing else."""
    import os
    place_dir = Path(place_dir)
    plan = json.load(open(place_dir / 'plan.json'))
    mine = chunk_jobs(plan['jobs'], chunk, of)
    log(f'chunk {chunk}/{of}: {len(mine)} shard placements')
    out = place_dir / 'joins' / f'{chunk:03d}-of-{of:03d}.jsonl'
    out.parent.mkdir(exist_ok=True)
    tmp = Path(f'{out}.tmp')
    with open(tmp, 'w') as f:
        for r in run_jobs(plan, mine, procs, regions_root, mask_cache=False, log=log):
            f.write(json.dumps(r, separators=(',', ':')) + '\n')
    os.replace(tmp, out)


def read_results(place_dir):
    """Every result of every chunk file under PLACE_DIR/joins (all chunks of one --of)."""
    files = sorted((Path(place_dir) / 'joins').glob('*.jsonl'))
    ofs = {f.name.split('-of-')[1] for f in files}
    if len(ofs) != 1:
        raise SystemExit(f'{place_dir}/joins: chunk files of different --of: {sorted(ofs)}')
    for f in files:
        for line in open(f):
            yield json.loads(line)


def reduce(plan, results, group_at, levels, log=print):
    """Placement, step 3 (in the release): pool every job's join levels per query (median over the
    shards holding both), assign each query (insert.assign), and return the rows of placements.tsv
    and the status counts. Every planned job must have exactly one result."""
    from . import insert
    want = {(q, shard) for q, _, shard in plan['jobs']}
    seen = set()
    pooled = collections.defaultdict(lambda: collections.defaultdict(list))
    shards, moved = collections.Counter(), collections.Counter()
    for r in results:
        key = (r['q'], r['shard'])
        if key not in want or key in seen:
            raise SystemExit(f'placement result {key} is {"not planned" if key not in want else "duplicated"}')
        seen.add(key)
        shards[r['q']] += 1
        moved[r['q']] += r['moved']
        for y, v in r['joins'].items():
            pooled[r['q']][y].append(v)
    if seen != want:
        raise SystemExit(f'{len(want - seen)} planned shard placements have no result (e.g. {sorted(want - seen)[:3]})')
    rows = [dict(r) for r in plan['rows']]
    for row in rows:
        q = row['input_id']
        if q not in plan['its2']:
            continue
        med = {y: float(np.median(v)) for y, v in pooled[q].items()}
        groups, avg = insert.assign(med, group_at, levels)
        k = next((i for i, g in enumerate(groups) if g is not None), None)
        row.update(status='placed' if k is not None else 'ungrouped',
                   group_id=groups[k] if k is not None else None, level=levels[k] if k is not None else None,
                   join=None if avg is None else round(avg, 5), shards=shards[q], rearranged=moved[q])
    status = collections.Counter(r['status'] for r in rows)
    log(f'ITS2-only placements: {dict(status)}')
    return rows, dict(status)


def check_plan(plan, comp_of, covers_of):
    """A placement plan fits this release: the same components, and the same ones built."""
    if plan['components_sha256'] != components_digest(comp_of):
        raise SystemExit('placement plan: made against other components than this release (another intake or --min-id)')
    if plan['built'] != sorted(covers_of):
        raise SystemExit(f'placement plan: built components differ from this release ({len(plan["built"])} vs {len(covers_of)})')


def place_inputs(intake_dir, classes, built_group_at, levels, comp_of, covers_of, threads=8, log=print):
    """Rows of placements.tsv for every ITS2-only input, all three steps in this process."""
    plan = plan_inputs(intake_dir, classes, comp_of, covers_of, threads, log)
    return reduce(plan, run_jobs(plan, plan['jobs'], threads, log=log), built_group_at, levels, log)
