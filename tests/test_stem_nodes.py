import numpy as np

from seqhesion import export
from seqhesion.coassoc import hierarchy as linkage
from seqhesion.hierarchy import LEVELS, linkage_nodes, stem_nodes

# four tips: 0 and 1 identical, 2 joins them at 0.003, 3 joins all at 0.02; three shards agree
PAIRS = {(0, 1): 0.0, (0, 2): 0.003, (1, 2): 0.003, (0, 3): 0.02, (1, 3): 0.02, (2, 3): 0.02}
N = 4


def votes():
    K = np.array(sorted(i * N + j for i, j in PAIRS for _ in range(3)))
    V = np.array([PAIRS[(k // N, k % N)] for k in K])
    return K, V


def tree():
    keys = np.array(sorted(i * N + j for i, j in PAIRS))
    med = np.array([PAIRS[(k // N, k % N)] for k in keys])
    Z = linkage(keys, med, N)[0]
    return [[int(r[0]), int(r[1]), float(r[2]), int(r[3])] for r in Z]


def test_a_node_living_between_levels_is_a_stem_node():
    Z = tree()
    assert [(s, h) for _, s, h, _ in linkage_nodes(Z, N, 0.1)] == [([0, 1], 0.0), ([0, 1, 2], 0.003), ([0, 1, 2, 3], 0.02)]
    st = stem_nodes(Z, N, votes(), LEVELS)
    # {0,1} lives 0 - 0.003 (no level inside, stem 0.003); {0,1,2} spans 0.005 (a level group)
    assert [(x['tips'], x['measured_at']) for x in st] == [([0, 1], 0.0025)]
    assert st[0]['cohesion'] == 1.0 and st[0]['held'] == 1.0


def test_a_short_stem_is_not_a_stem_node():
    assert stem_nodes(tree(), N, votes(), LEVELS, stem_min=0.004) == []


def test_parent_is_the_nearest_enclosing_group_or_node():
    Z = tree()
    tips = ['a', 'b', 'c', 'd']
    rows = [{'group_id': 'L1', 'kind': 'level', 'parent_id': 'L2'},       # {a,b,c}, level 0.005
            {'group_id': 'L2', 'kind': 'level', 'parent_id': None},       # all, level 0.02
            {'group_id': 'N1', 'kind': 'node', 'parent_id': None}]        # {a,b}, a stem node
    tips_of = {'L1': ['a', 'b', 'c'], 'L2': tips, 'N1': ['a', 'b']}
    export.node_columns(rows, tips_of, [({'Z': Z}, tips)])
    by = {r['group_id']: r for r in rows}
    assert by['N1']['parent_id'] == 'L1' and by['L1']['parent_id'] == 'L2' and by['L2']['parent_id'] is None
    assert by['N1']['dendrogram_node'] == 'm0' and by['N1']['height_top'] == 0.003 and by['N1']['linkage_stem'] == 0.003
