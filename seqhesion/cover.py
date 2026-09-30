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


def plan(nbrs, all_ids, scaffold, covers=2, depth=3, size=150, quota=0, quota_id=97.0, ident=None, log=print):
    """The manifest of `covers` independent covers. With a context quota, every shard gets at
    least `quota` members below `quota_id` identity to its seed (its nearest such neighbours
    appended), since a dense species otherwise fills a window and is never seen against outsiders."""
    import random
    manifest = {'scaffold': scaffold, 'shards': {}, 'quota': {'n': quota, 'below_identity': quota_id} if quota else None}
    topped = 0
    for ci in range(covers):
        rng = random.Random(100 + ci)
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
    """A shard's tree input, in order: its members sorted, then the scaffold tips not already in it."""
    members = manifest['shards'][name]['members']
    mem = set(members)
    return sorted(members) + [s for s in manifest['scaffold'] if s not in mem]


_SEQS = {}


def _tree_job(job):
    from . import tools
    name, ids, tree_dir = job
    return name, tools.cached_build_tree({i: _SEQS[i] for i in ids}, tree_dir, name, order=ids)[1]


def _init_seqs(seqs):
    _SEQS.update(seqs)


def build_trees(manifest, seqs, tree_dir, procs=6, log=print):
    """One tree per shard under tree_dir, reusing any whose key matches. Returns trees reused."""
    from multiprocessing import Pool
    jobs = [(name, shard_input(manifest, name), str(tree_dir)) for name in manifest['shards']]
    reused = 0
    with Pool(procs, initializer=_init_seqs, initargs=(seqs,)) as pool:
        for i, (_, r) in enumerate(pool.imap_unordered(_tree_job, jobs), 1):
            reused += r
            if i % 40 == 0:
                log(f'  {i}/{len(jobs)} trees')
    log(f'trees: {reused} reused, {len(jobs) - reused} built')
    return reused
