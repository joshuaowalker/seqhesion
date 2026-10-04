"""Invariants of the cover builder: no repeated shards, deterministic, every input accounted for."""
import random

from seqhesion import cover


def ring(n, k):
    """n sequences on a ring; each one's neighbours are the k nearest on the ring."""
    ids = [f's{i:03d}' for i in range(n)]
    nbrs = {}
    for i, s in enumerate(ids):
        order = sorted((min((i - j) % n, (j - i) % n), j) for j in range(n) if j != i)
        nbrs[s] = [ids[j] for _, j in order[:k]]
    return ids, nbrs


def test_no_repeated_seed_or_member_set():
    ids, nbrs = ring(120, 30)
    shards = cover.greedy_cover(nbrs, 20, 3, random.Random(1))
    assert len({s for s, _ in shards}) == len(shards)
    assert len({frozenset(m) for _, m in shards}) == len(shards)
    assert all(len(m) == 20 and len(set(m)) == 20 for _, m in shards)


def test_deterministic_for_a_given_seed():
    ids, nbrs = ring(120, 30)
    assert cover.greedy_cover(nbrs, 20, 3, random.Random(7)) == cover.greedy_cover(nbrs, 20, 3, random.Random(7))


def test_every_input_is_accounted_for():
    ids, nbrs = ring(120, 30)
    nbrs['lonely'] = []                              # no neighbours at all
    nbrs['sparse'] = ids[:6]                         # too few to seed a full shard
    everyone = set(ids) | {'lonely', 'sparse'}
    rng = random.Random(3)
    shards = cover.greedy_cover(nbrs, 20, 3, rng)
    acct = cover.account(everyone, nbrs, shards, 20, 3)
    assert acct['covered'] + len(acct['under']) + len(acct['isolated']) == len(everyone)
    assert 'lonely' in acct['isolated'] and acct['under']['sparse'] == 0
    shards += cover.sparse_cover(nbrs, everyone, shards, 20, 3, rng)
    acct = cover.account(everyone, nbrs, shards, 20, 3)
    assert acct['under'].get('sparse', 3) >= 1       # now in at least one tree
    assert 'lonely' in acct['isolated']              # still an explicit outcome, not an omission


def test_identical_neighbourhoods_are_reported_not_hidden():
    """Five sequences that are each other's only neighbours plus a big ring: the five can
    only ever form one distinct shard, so depth 3 is unreachable and must be reported."""
    ids, nbrs = ring(60, 25)
    clique = [f'q{i}' for i in range(5)]
    for q in clique:
        nbrs[q] = [x for x in clique if x != q]
    everyone = set(ids) | set(clique)
    rng = random.Random(5)
    shards = cover.greedy_cover(nbrs, 20, 3, rng)
    shards += cover.sparse_cover(nbrs, everyone, shards, 20, 3, rng)
    assert len({frozenset(m) for _, m in shards}) == len(shards)
    acct = cover.account(everyone, nbrs, shards, 20, 3)
    assert all(acct['under'][q] == 1 for q in clique)


def test_as_context_and_shard_input():
    from seqhesion.cover import as_context, shard_input
    man = {'scaffold': ['s1', 'b', 's2'], 'shards': {
        'c0.0000': {'cover': 0, 'seed': 'a', 'members': ['a', 'x', 'b', 'y']},
        'c0.0001': {'cover': 0, 'seed': 'a', 'members': ['a', 'x', 'y']}}}
    small = as_context(man, {'a', 'b'})
    assert list(small['shards']) == ['c0.0000']                      # one own tip left: dropped
    s = small['shards']['c0.0000']
    assert s['members'] == ['a', 'b'] and s['context'] == ['x', 'y']
    # the tree input is the same as before the split: everything sorted, then the scaffold not in it
    assert shard_input(small, 'c0.0000') == shard_input(man, 'c0.0000') == ['a', 'b', 'x', 'y', 's1', 's2']


def test_twins_are_cover1_shards_with_a_cover0_member_set():
    man = {'scaffold': ['s'], 'shards': {
        'a0': {'cover': 0, 'members': ['x', 'y', 'z']},
        'b0': {'cover': 0, 'members': ['p', 'q'], 'context': ['c']},
        'a1': {'cover': 1, 'members': ['z', 'x', 'y']},              # same set, other order: a twin
        'b1': {'cover': 1, 'members': ['p', 'q']},                   # context differs: not a twin
        'c1': {'cover': 1, 'members': ['x', 'y']}}}                  # a subset: not a twin
    assert cover.twins(man) == {'a1': 'a0'}
