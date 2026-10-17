"""Overlapping shard covers.

A shard is a seed plus its nearest neighbours. A cover is a set of shards such
that every eligible sequence lies in at least `depth` of them. Seeds are chosen
greedily from the currently least-covered sequences; `rng` only breaks ties, so
different rng seeds give different but equally valid covers. No seed and no
member set is used twice: a repeated shard counts as a replicate without being
one (the first covers of one test region had 102 such repeats among 311 shards).

Sequences with fewer than size - 1 neighbours above the search threshold cannot
seed a full shard. greedy_cover() alone leaves them out (4% of that region); sparse_cover()
then gives each of them smaller shards made of whatever neighbours they have.
account() reconciles every input sequence against a cover, so that nothing is
silently missing from the metrics: each sequence is 'covered' (in >= depth shards),
'under' (in fewer; listed with its count) or 'isolated' (too few neighbours to be
in any tree at all; still an outcome — a singleton nobody has placed).
"""
import collections
import csv
import json
import os
from pathlib import Path


def load_knn(path):
    """{query: [targets, most similar first]} from a vsearch userout (query, target, id)."""
    nbrs = collections.defaultdict(list)
    with open(path) as f:
        for q, t, ident in csv.reader(f, delimiter='\t'):
            nbrs[q].append((float(ident), t))
    return {q: [t for _, t in sorted(v, reverse=True)] for q, v in nbrs.items()}


def shard(seed, nbrs, size):
    return [seed] + nbrs[seed][:size - 1]


def greedy_cover(nbrs, size, depth, rng):
    """Returns [(seed, members)], covering every sequence that has a full
    neighbourhood at least `depth` times."""
    eligible = sorted(s for s, v in nbrs.items() if len(v) >= size - 1)
    coverage = dict.fromkeys(eligible, 0)
    shards, used, seen = [], set(), set()
    while True:
        low = min(coverage.values(), default=depth)
        if low >= depth:
            return shards
        cands = [s for s in eligible if coverage.get(s) == low and s not in used]
        if cands:
            seed = rng.choice(cands)
            members = shard(seed, nbrs, size)
        else:
            # Every least-covered sequence has already seeded a shard, so it is in nobody's
            # neighbourhood but its own. Seeding it again would repeat the same shard (a
            # replicate that is not one); instead seed its nearest unused neighbour that
            # gives a new shard and put the sequence into it, so it is seen in a different
            # context.
            x = rng.choice([s for s in eligible if coverage.get(s) == low])
            for seed in nbrs[x]:
                if seed in used or seed not in coverage:
                    continue
                members = shard(seed, nbrs, size)
                if x not in members:
                    members = members[:-1] + [x]
                if frozenset(members) not in seen:
                    break
            else:                          # no fresh context left: leave it under-covered
                del coverage[x]
                continue
        used.add(seed)
        if frozenset(members) in seen:     # another seed with exactly this neighbourhood
            continue
        seen.add(frozenset(members))
        shards.append((seed, members))
        for m in members:
            if m in coverage:
                coverage[m] += 1


def sparse_cover(nbrs, all_ids, shards, size, depth, rng, min_size=4):
    """Extra shards for the sequences greedy_cover() cannot seed: seed + all of its
    neighbours (fewer than size - 1), provided that makes at least min_size sequences.
    Run after greedy_cover so the full-size shards are unaffected. Returns the new
    [(seed, members)] only."""
    coverage = collections.Counter(m for _, m in shards for m in m)
    seen = {frozenset(m) for _, m in shards}
    sparse = sorted(s for s in all_ids if min_size - 1 <= len(nbrs.get(s, ())) < size - 1)
    new, used = [], set()
    while True:
        low = [s for s in sparse if coverage[s] < depth and s not in used]
        if not low:
            return new
        least = min(coverage[s] for s in low)
        seed = rng.choice([s for s in low if coverage[s] == least])
        used.add(seed)
        members = shard(seed, nbrs, size)
        if frozenset(members) in seen:
            continue
        seen.add(frozenset(members))
        new.append((seed, members))
        coverage.update(members)


def account(all_ids, nbrs, shards, size, depth, min_size=4):
    """Reconcile every input sequence against a cover (see module docstring)."""
    coverage = collections.Counter(m for _, m in shards for m in m)
    out = {'input': len(all_ids), 'shards': len(shards), 'depth': depth,
           'covered': 0, 'under': {}, 'isolated': {}}
    for s in sorted(all_ids):
        if coverage[s] >= depth:
            out['covered'] += 1
        elif coverage[s] == 0 and len(nbrs.get(s, ())) < min_size - 1:
            out['isolated'][s] = len(nbrs.get(s, ()))
        else:
            out['under'][s] = coverage[s]
    assert out['covered'] + len(out['under']) + len(out['isolated']) == len(all_ids)
    return out


# ---------------------------------------------------------------------------------------------
# Building a region's covers: shards (seed + nearest neighbours + a shared scaffold), one tree each.

N_SCAFFOLD = 30


def load_identity(path):
    """{query: {target: identity}} from the same vsearch userout as load_knn."""
    ident = {}
    for line in open(path):
        q, t, i = line.rstrip('\n').split('\t')
        ident.setdefault(q, {})[t] = float(i)
    return ident


def scaffold_of(centroids_fasta, n=N_SCAFFOLD):
    """The n most abundant 85% centroids (headers `id;size=N`): rooting context shared by every
    shard, pruned again before any co-association is counted."""
    import re
    from .fasta import read_fasta
    sized = [(int(re.search(r'size=(\d+)', h).group(1)), re.sub(r';size=\d+;?$', '', h)) for h in read_fasta(centroids_fasta)]
    return [c for _, c in sorted(sized, reverse=True)[:n]]


def plan(nbrs, all_ids, scaffold, covers=2, depth=3, size=150, quota=0, quota_id=97.0, ident=None, log=print, draw=0):
    """The manifest of `covers` independent covers. With a context quota, every shard gets at
    least `quota` members below `quota_id` identity to its seed (its nearest such neighbours
    appended), since a dense species otherwise fills a window and is never seen against outsiders."""
    import random
    manifest = {'scaffold': scaffold, 'shards': {}, 'quota': {'n': quota, 'below_identity': quota_id} if quota else None}
    topped = 0
    for ci in range(covers):
        rng = random.Random(100 + ci + 1000 * draw)
        shards = greedy_cover(nbrs, size, depth, rng)
        n_full = len(shards)
        shards += sparse_cover(nbrs, all_ids, shards, size, depth, rng)
        acct = account(all_ids, nbrs, shards, size, depth)
        manifest.setdefault('coverage', {})[str(ci)] = acct
        log(f'cover {ci}: {n_full} full shards + {len(shards) - n_full} for sparse sequences; of {acct["input"]} input sequences '
            f'{acct["covered"]} are in >= {depth} shards, {len(acct["under"])} in fewer, {len(acct["isolated"])} isolated')
        for k, (seed, members) in enumerate(shards):
            if quota:
                far = [t for t in members if t != seed and ident.get(seed, {}).get(t, 0.0) < quota_id]
                if len(far) < quota:
                    have = set(members)
                    extra = [t for t in nbrs.get(seed, []) if t not in have and ident[seed].get(t, 0.0) < quota_id][:quota - len(far)]
                    if extra:
                        members = list(members) + extra
                        topped += 1
            manifest['shards'][f'c{ci}.{k:04d}'] = {'cover': ci, 'seed': seed, 'members': list(members)}
    if quota:
        log(f'{topped} shards topped up to the context quota')
    return manifest


def shard_input(manifest, name):
    """A shard's tree input, in order: its members and context tips sorted together, then the
    scaffold tips not already in it. (Only small regions' shards have context; see as_context.)"""
    s = manifest['shards'][name]
    inside = list(s['members']) + list(s.get('context', []))
    have = set(inside)
    return sorted(inside) + [t for t in manifest['scaffold'] if t not in have]


def twins(manifest):
    """{cover-1 shard: the cover-0 shard with exactly its member set}. The two covers are drawn
    independently, so in dense neighbourhoods they often draw the same shard (55-61% of cover 1
    over the Agaricales build): the same sample, so the same tree. The hierarchy pools both covers
    and counts such a shard once (Josh, 2026-10-04); placement skips it likewise."""
    key = lambda s: frozenset(s['members']) | frozenset(s.get('context', []))  # noqa: E731
    first = {key(s): name for name, s in sorted(manifest['shards'].items()) if s['cover'] == 0}
    return {name: first[key(s)] for name, s in sorted(manifest['shards'].items())
            if s['cover'] == 1 and key(s) in first}


def as_context(manifest, own):
    """A small region's manifest: each shard's tips from outside the component (`own`) become
    context -- in the tree, rooting candidates, pruned before co-association -- and shards left with
    fewer than two of the component's own tips are dropped (no pair to join)."""
    shards = {}
    for name, s in manifest['shards'].items():
        mine = [m for m in s['members'] if m in own]
        if len(mine) >= 2:
            shards[name] = {**s, 'members': mine, 'context': [m for m in s['members'] if m not in own]}
    return {**manifest, 'shards': shards, 'small': True}


_SEQS = {}


def _tree_job(job):
    from . import tools
    name, ids, tree_dir = job
    return name, tools.cached_build_tree({i: _SEQS[i] for i in ids}, tree_dir, name, order=ids)[1]


def _init_seqs(seqs):
    _SEQS.update(seqs)


TREE_FILES = ('.fasta', '.aln.fasta', '.trim.fasta', '.nwk', '.nwk.key.json')     # the key sidecar last


def plan_region(region_dir, depth=3, size=150, quota=30, quota_id=97.0, covers=2, log=print):
    """The covers manifest of a region dir and the sequences its trees need ({id: seq}): small regions
    draw on their corpus-wide context (seeds from their own tips only)."""
    from .fasta import read_fasta
    rd = Path(region_dir)
    seqs = read_fasta(rd / 'comp.fasta')
    own = set(seqs)
    small = json.load(open(rd / 'region.json')).get('kind') == 'small'
    nbrs = load_knn(rd / 'knn.tsv')
    if small:
        seqs.update(read_fasta(rd / 'context.fasta'))
        nbrs = {t: [x for x in nbrs.get(t, []) if x != t] for t in sorted(own)}
    ident = load_identity(rd / 'knn.tsv') if quota else None
    man = plan(nbrs, own, scaffold_of(rd / 'comp_cent_0.85.fasta'), covers, depth, size, quota, quota_id, ident, log=log)
    if small:
        man = as_context(man, own)
    return man, seqs


def region_seqs(region_dir):
    """Every sequence a region's shard trees can hold: its own tips and, for a small region, its context."""
    from .fasta import read_fasta
    rd = Path(region_dir)
    seqs = read_fasta(rd / 'comp.fasta')
    if (rd / 'context.fasta').exists():
        seqs.update(read_fasta(rd / 'context.fasta'))
    return seqs


def distinct_inputs(manifest):
    """{tree input (ids in order): [shard names]}: shards with the same input are one tree (twins)."""
    out = {}
    for name in sorted(manifest['shards']):
        out.setdefault(tuple(shard_input(manifest, name)), []).append(name)
    return out


def reuse_index(dirs):
    """{key: path stem} of every complete tree (all TREE_FILES present) under the given directories."""
    idx = {}
    for d in map(Path, dirs):
        if not d.is_dir():
            continue
        for f in sorted(d.glob('*.nwk.key.json')):
            stem = d / f.name[:-len('.nwk.key.json')]
            if all(Path(f'{stem}{ext}').exists() for ext in TREE_FILES):
                idx.setdefault(json.load(open(f))['key'], stem)
    return idx


def link_reusable(manifest, seqs, tree_dir, reuse_dirs, log=print):
    """Hard-link into tree_dir, under the shard's own name, every earlier tree whose key is exactly this
    shard's (the same sequences in the same order, recipe and tools: seqhesion.cache), wherever it was
    built and whatever it was called: a cover redrawn after a few sequences were added keeps many shards
    under other names. Tree files are never modified in place, so a link cannot change the source.
    Returns how many trees were linked."""
    from . import tools
    idx = reuse_index(reuse_dirs)
    td = Path(tree_dir)
    td.mkdir(parents=True, exist_ok=True)
    linked = 0
    for ids, names in distinct_inputs(manifest).items():
        name = names[0]
        if Path(f'{td / name}.nwk.key.json').exists():
            continue                                     # already here: cache.valid decides when building
        src = idx.get(tools.build_key({i: seqs[i] for i in ids}, order=list(ids)))
        if src is None:
            continue
        for ext in TREE_FILES:                           # the key sidecar last
            tmp = Path(f'{td / name}{ext}.tmp')
            if tmp.exists():
                tmp.unlink()
            os.link(f'{src}{ext}', tmp)
            os.replace(tmp, f'{td / name}{ext}')
        linked += 1
    log(f'reuse: {linked} trees linked from {len(idx)} earlier ones')
    return linked


def build_trees(manifest, seqs, tree_dir, procs=6, log=print, reuse=()):
    """One tree per shard under tree_dir, reusing any whose key matches. Returns trees reused.

    Shards with the same input (ids in the same order; ids are content hashes, so the same
    sequences) are built once: in dense neighbourhoods the two covers often draw identical shards
    (30% of v4's trees). The twin gets byte copies of the first one's files, which is exactly
    what building it again would give (one invocation, one result: verified under --thread 1)."""
    import shutil
    from multiprocessing import Pool
    from . import cache
    if reuse:
        link_reusable(manifest, seqs, tree_dir, reuse, log=log)
    twins = distinct_inputs(manifest)
    jobs = [(names[0], list(ids), str(tree_dir)) for ids, names in twins.items()]
    reused = 0
    with Pool(procs, initializer=_init_seqs, initargs=(seqs,)) as pool:
        for i, (_, r) in enumerate(pool.imap_unordered(_tree_job, jobs), 1):
            reused += r
            if i % 40 == 0:
                log(f'  {i}/{len(jobs)} trees')
    td = Path(tree_dir)
    copied = 0
    for ids, names in twins.items():
        first, key = names[0], json.load(open(td / f'{names[0]}.nwk.key.json'))['key']
        for twin in names[1:]:
            if cache.valid(td / f'{twin}.nwk', key, tips=ids):
                reused += 1
                continue
            for ext in TREE_FILES:
                shutil.copyfile(td / f'{first}{ext}', td / f'{twin}{ext}.tmp')
                os.replace(td / f'{twin}{ext}.tmp', td / f'{twin}{ext}')
            copied += 1
    n = len(manifest['shards'])
    log(f'trees: {reused} reused, {n - reused - copied} built, {copied} copied from an identical shard')
    return reused


def missing_trees(region_dirs, covers='covers', log=print):
    """The trees the planned covers of these regions still need (not present with the right key): a
    job list for build_jobs, one job per distinct tree input with every shard name it serves."""
    from . import cache, tools
    jobs = []
    for rd in map(Path, region_dirs):
        man = json.load(open(rd / covers / 'manifest.json'))
        seqs = region_seqs(rd)
        td = rd / covers / 'trees'
        for ids, names in distinct_inputs(man).items():
            key = tools.build_key({i: seqs[i] for i in ids}, order=list(ids))
            if not cache.valid(td / f'{names[0]}.nwk', key, tips=list(ids)):
                jobs.append({'region': rd.name, 'names': names, 'ids': list(ids)})
    log(f'{len(jobs)} trees to build over {len(region_dirs)} regions')
    return jobs


def chunk_jobs(jobs, chunk, of):
    """Chunk `chunk` of `of`: balanced by size squared (L-INS-i cost), largest first; deterministic."""
    order = sorted(range(len(jobs)), key=lambda k: (-len(jobs[k]['ids']), jobs[k]['region'], jobs[k]['names'][0]))
    load, mine = [0] * of, []
    for k in order:
        c = min(range(of), key=lambda i: (load[i], i))
        load[c] += len(jobs[k]['ids']) ** 2
        if c == chunk:
            mine.append(jobs[k])
    return mine


def build_jobs(jobs, regions_root, covers='covers', chunk=0, of=1, procs=6, log=print):
    """Build one chunk of a job list (missing_trees) under REGIONS_ROOT/<region>/<covers>/trees, then
    give each tree's twins (other shard names with the same input) byte copies, as build_trees does."""
    import shutil
    from multiprocessing import Pool
    mine = chunk_jobs(jobs, chunk, of)
    root = Path(regions_root)
    seqs = {}
    for r in sorted({j['region'] for j in mine}):
        seqs.update(region_seqs(root / r))
    log(f'chunk {chunk}/{of}: {len(mine)} trees')
    work = [(j['names'][0], j['ids'], str(root / j['region'] / covers / 'trees')) for j in mine]
    for w in work:
        Path(w[2]).mkdir(parents=True, exist_ok=True)
    with Pool(procs, initializer=_init_seqs, initargs=(seqs,)) as pool:
        for i, _ in enumerate(pool.imap_unordered(_tree_job, work), 1):
            if i % 40 == 0:
                log(f'  {i}/{len(work)} trees')
    copied = 0
    for j in mine:
        td = root / j['region'] / covers / 'trees'
        for twin in j['names'][1:]:
            for ext in TREE_FILES:
                shutil.copyfile(td / f"{j['names'][0]}{ext}", td / f'{twin}{ext}.tmp')
                os.replace(td / f'{twin}{ext}.tmp', td / f'{twin}{ext}')
            copied += 1
    log(f'chunk {chunk}/{of}: {len(mine)} trees built, {copied} twins copied')
