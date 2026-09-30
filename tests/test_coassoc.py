import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from seqhesion.coassoc import UNJOINED, forest_roots, observed_average_linkage, pooled, separation, shard_evidence, votes_at


def pairs(D, keep=None):
    n = len(D)
    i, j = np.triu_indices(n, 1)
    m = np.ones(len(i), bool) if keep is None else keep[i, j]
    return i[m], j[m], D[i, j][m]


def same_partition(a, b):
    return len({(x, y) for x, y in zip(a, b)}) == len(set(a)) == len(set(b))


def test_all_pairs_observed_is_ordinary_average_linkage():
    rng = np.random.default_rng(1)
    X = rng.random((40, 3))
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    Z, obs, pos = observed_average_linkage(*pairs(D), 40)
    ref = linkage(squareform(D, checks=False), method='average')
    assert np.allclose(Z[:, 2], ref[:, 2])
    assert np.array_equal(obs, pos)                 # every cross pair observed
    for t in (0.2, 0.4, 0.8):
        assert same_partition(fcluster(Z, t, 'distance'), fcluster(ref, t, 'distance'))


def test_unobserved_pair_is_no_evidence_not_distance():
    # 0-1 close, 1-2 at 0.2, 0-2 never observed: 2 joins at 0.2 (mean over the one observed
    # cross pair); filling the missing pair with 1.0 would put it at (0.2 + 1.0) / 2 = 0.6
    Z, obs, pos = observed_average_linkage([0, 1], [1, 2], [0.1, 0.2], 3)
    assert Z[:, 2].tolist() == [0.1, 0.2]
    assert obs.tolist() == [1, 1] and pos.tolist() == [1, 2]


def test_no_cross_evidence_leaves_a_forest():
    Z, obs, _ = observed_average_linkage([0, 2, 2], [1, 3, 4], [0.1, 0.1, 0.3], 5)
    assert forest_roots(Z) == 2 and Z[-1, 2] == UNJOINED and obs[-1] == 0
    assert Z[-1, 3] == 5
    p = fcluster(Z, 10.0, 'distance')
    assert len(set(p)) == 2 and p[0] == p[1] != p[2] == p[3] == p[4]


def test_sparse_heights_monotone_and_order_free():
    rng = np.random.default_rng(2)
    n = 120
    X = rng.random((n, 2))
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    keep = D < 0.25                                   # neighbourhoods, like shards
    i, j, d = pairs(D, keep)
    Z, _, _ = observed_average_linkage(i, j, d, n)
    h = Z[Z[:, 2] != UNJOINED, 2]
    assert np.all(np.diff(h) >= 0)
    o = rng.permutation(len(i))
    Z2, _, _ = observed_average_linkage(i[o], j[o], d[o], n)
    assert np.array_equal(Z, Z2)


def test_votes_at_counts_shards_per_pair():
    K = np.array([1, 1, 1, 5, 7, 7])
    V = np.array([0.01, 0.02, 0.09, 0.2, 0.01, 0.01])
    k, h, v = votes_at(K, V, 0.05)
    assert k.tolist() == [1, 5, 7] and h.tolist() == [2, 0, 2] and v.tolist() == [3, 1, 2]


def test_pooled_weighs_votes_and_pulls_by_group():
    # tips 0,1,2 in group 0; 3,4 in group 1; 5 ungrouped
    n = 6
    g = np.array([0, 0, 0, 1, 1, -1])
    key = lambda a, b: min(a, b) * n + max(a, b)
    pairs = {(0, 1): (10, 10), (0, 2): (1, 10), (1, 2): (10, 10),     # (hits, votes)
             (0, 3): (1, 1),                                          # one shard, joined: raw 1.00
             (0, 4): (0, 9),
             (2, 5): (0, 4), (3, 4): (5, 5)}
    ks = sorted(pairs, key=lambda p: key(*p))
    keys = np.array([key(*p) for p in ks])
    hits = np.array([pairs[p][0] for p in ks])
    votes = np.array([pairs[p][1] for p in ks])
    s = pooled(keys, hits, votes, n, g)
    assert s['membership'][0] == 11 / 20 and s['m_votes'][0] == 20 and s['m_partners'][0] == 2
    assert s['cohesion'][0] == 21 / 30 and s['c_votes'][0] == 30 and s['c_pairs'][0] == 3
    # tip 0's single-shard 1.00 with tip 3 is pooled with the 9 votes against tip 4: 1 / 10
    assert s['pull'][0] == 0.1 and s['p_group'][0] == 1 and s['p_votes'][0] == 10 and s['p_partners'][0] == 2
    # tip 2 sees outside only tip 5, never joined: none at this level, with its evidence
    assert s['p_group'][2] == -3 and s['p_votes'][2] == 4
    # tip 1 has no outside pair at all
    assert s['p_group'][1] == -2
    # tip 5 (ungrouped) pulled towards group 0 with nothing joining
    assert s['p_group'][5] == -3 and np.isnan(s['membership'][5])


def test_shard_evidence_verdicts_stems_and_hits():
    from ete3 import Tree
    index = {k: i for i, k in enumerate('abcde')}
    t = Tree('(((a:0.01,b:0.01):0.02,c:0.03):0.1,(d:0.01,e:0.01):0.1);', format=5)
    m = lambda s: sum(1 << index[k] for k in s)  # noqa: E731
    groups = [m('ab'), m('abc'), m('bc'), m('ad'), m('abcde')]
    of = [[g for g, gm in enumerate(groups) if gm >> i & 1] for i in range(5)]
    rows = {r[0]: r for r in shard_evidence(t, index, 500, of, groups, [0] * 5, [0.03])}
    assert rows[0][3] == 'clade' and abs(rows[0][4] - 10.0) < 1e-9 and rows[0][5:] == (1, 1)
    assert rows[1][3] == 'clade' and rows[1][5:] == (1, 3)      # only a-b joined within 0.03
    assert rows[2][3] == 'conflict'                              # (a,b) takes b with outsider a
    assert rows[3][3] == 'conflict'
    assert rows[4][3] == 'whole' and rows[4][4] is None


def test_separation_names_nearest_outside_groups_not_subgroups():
    n = 7
    key = lambda a, b: min(a, b) * n + max(a, b)  # noqa: E731
    d = {(0, 1): 0.01, (0, 2): 0.02, (1, 2): 0.015, (2, 3): 0.05, (0, 4): 0.08, (1, 5): 0.03,
         (2, 6): 0.06}
    ks = sorted(d, key=lambda p: key(*p))
    keys = np.array([key(*p) for p in ks])
    dist = np.array([d[p] for p in ks])
    group_of = np.array([0, 0, 0, 1, 1, -1, 2])  # group {0,1,2}; relatives: groups 1 and 2, stray 5
    s = separation(keys, dist, n, [0, 1, 2], group_of)
    assert s['spread'] == 0.02 and s['within_pairs'] == 3
    assert s['nearest'] == (0.05, 3, 1)          # group 1 via tip 3 (not tip 4 at 0.08)
    assert s['second'] == (0.06, 6, 2)
    assert s['stray'] == (0.03, 5) and s['strays_closer'] == 1


def test_compatible_join_reads_a_polytomy_as_missing_information():
    from ete3 import Tree
    from seqhesion.coassoc import join_matrix
    index = {k: i for i, k in enumerate('abcd')}
    # a and b identical, c and d far, all four in one polytomy
    t = Tree('(a:0.0,b:0.0,c:0.02,d:0.03);', format=5)
    ids, J, _, _ = join_matrix(t, index, 'diameter')
    k = {x: ids.index(index[x]) for x in 'abcd'}
    assert J[k['a'], k['b']] == 0.05                      # the polytomy's full width
    ids, J, _, _ = join_matrix(t, index, 'compatible')
    assert J[k['a'], k['b']] == 0.0 and J[k['a'], k['c']] == 0.02 and J[k['c'], k['d']] == 0.05
    # a resolved node is the same either way
    t = Tree('((a:0.01,b:0.01):0.02,(c:0.01,d:0.01):0.02);', format=5)
    assert np.array_equal(join_matrix(t, index, 'diameter')[1], join_matrix(t, index, 'compatible')[1])


def test_identical_join_only_joins_tips_below_resolution_inside_a_polytomy():
    from ete3 import Tree
    from seqhesion.coassoc import join_matrix
    index = {k: i for i, k in enumerate('abcd')}
    t = Tree('(a:0.0,b:0.001,c:0.02,d:0.03);', format=5)
    ids, J, _, _ = join_matrix(t, index, 'identical', floor=0.002)
    k = {x: ids.index(index[x]) for x in 'abcd'}
    assert J[k['a'], k['b']] == 0.001                     # indistinguishable: joined at their span
    assert J[k['a'], k['c']] == 0.05                      # distinguishable: the polytomy's width
