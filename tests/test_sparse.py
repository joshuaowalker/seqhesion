from scipy.cluster.hierarchy import fcluster
import numpy as np

from seqhesion import sparse


def _cut(Z, L):
    f = fcluster(np.array(Z, float), t=L, criterion='distance')
    g = {}
    for i, x in enumerate(f):
        g.setdefault(x, set()).add(i)
    return sorted(sorted(s) for s in g.values() if len(s) >= 2)


def test_combined_linkage_keeps_units_and_raises_low_layer_merges():
    # tips 0,1 form a unit at 0.1; 2 and 3 are lone tips. The fine linkage chains them above 0.1;
    # the layer joins units 1 and 2 at 0.08 (below the unit level) and then everything at 0.3
    fine = {'Z': [[0, 1, 0.05, 2], [4, 2, 0.2, 3], [5, 3, 0.25, 4]], 'observed': [1, 2, 3], 'possible': [1, 2, 3]}
    layer = {'units': 3, 'unit_of_tip': [0, 0, 1, 2], 'linkage': [[1, 2, 0.08, 2], [3, 0, 0.3, 3]],
             'observed': [1, 2], 'possible': [1, 2]}
    Z = sparse.combined_linkage(fine, layer, 4)['Z']
    assert len(Z) == 3
    assert _cut(Z, 0.1) == [[0, 1]]                 # the unit level is exactly the fine one
    assert _cut(Z, 0.15) == [[0, 1], [2, 3]]        # the layer's low merge sits just above it
    assert _cut(Z, 0.4) == [[0, 1, 2, 3]]
