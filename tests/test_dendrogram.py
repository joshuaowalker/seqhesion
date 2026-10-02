import numpy as np
from ete3 import Tree

from seqhesion import coassoc, export


def linkage_of(pairs, n):
    i, j, d = (np.array(x) for x in zip(*pairs))
    Z, obs, pos = coassoc.observed_average_linkage(i, j, d.astype(float), n)
    return {'Z': [[int(r[0]), int(r[1]), None if not np.isfinite(r[2]) else float(r[2]), int(r[3])] for r in Z],
            'observed': [int(x) for x in obs], 'possible': [int(x) for x in pos]}


def test_one_tree_is_ultrametric_with_every_tip():
    tips = ['a', 'b', 'c', 'd']
    lk = linkage_of([(0, 1, 0.01), (2, 3, 0.02), (0, 2, 0.08), (1, 3, 0.1), (0, 3, 0.09), (1, 2, 0.07)], 4)
    lines, rows = export.dendrogram(lk, tips)
    assert len(lines) == 1 and len(rows) == 3
    t = Tree(lines[0], format=1)
    assert sorted(t.get_leaf_names()) == tips
    root_h = max(r['height'] for r in rows)
    for leaf in t.iter_leaves():
        assert abs(t.get_distance(leaf) - root_h) < 1e-9          # ultrametric: every tip at the root height
    assert abs(t.get_distance('a', 'b') - 2 * 0.01) < 1e-9


def test_a_forest_gives_one_tree_per_part_largest_first():
    tips = ['a', 'b', 'c', 'x', 'y']
    lk = linkage_of([(0, 1, 0.01), (1, 2, 0.02), (0, 2, 0.03), (3, 4, 0.005)], 5)   # no pair links {a,b,c} to {x,y}
    lines, rows = export.dendrogram(lk, tips)
    assert len(lines) == 2
    assert sorted(Tree(lines[0], format=1).get_leaf_names()) == ['a', 'b', 'c']
    assert sorted(Tree(lines[1], format=1).get_leaf_names()) == ['x', 'y']
    assert any(r['height'] is None for r in rows)


def test_a_singleton_part_is_its_own_line():
    tips = ['a', 'b', 'z']
    lk = linkage_of([(0, 1, 0.01)], 3)
    lines, _ = export.dendrogram(lk, tips)
    assert sorted(lines) == sorted(['(a:0.01,b:0.01)m0;', 'z:0;'])
