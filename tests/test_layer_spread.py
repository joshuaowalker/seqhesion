"""Spread of coarse groups (schema 0.8, 2026-10-18): measured as for fine groups, from the
component's local trees (patristic distance, median over the trees holding a pair). A broad group
whose few far members set its widest pair shows a small spread_median beside a large spread."""
import itertools
import json
import statistics

from ete3 import Tree

from seqhesion import sparse

# unit 0 = a1 ... a4; unit 1 = b1, b2, b3, with b3 alone far out on a long branch
NWK = ['((a1:0.01,a2:0.012,a3:0.011,a4:0.01):0.002,(b1:0.01,b2:0.012,b3:0.06):0.002);',
       '((a1:0.012,a2:0.01,a3:0.01,a4:0.011):0.002,(b1:0.012,b2:0.01,b3:0.05):0.002);']
TIPS = ['a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'b3']


def covers(tmp_path):
    cd = tmp_path / 'covers'
    (cd / 'coassoc').mkdir(parents=True)
    (cd / 'trees').mkdir()
    shards = {f'c{c}.0000': {'cover': c, 'seed': 'a1', 'members': TIPS} for c in (0, 1)}
    shards['c1.0000']['context'] = ['x']              # not a twin of c0.0000: both trees count
    json.dump({'scaffold': [], 'shards': shards}, open(cd / 'manifest.json', 'w'))
    for c, nwk in enumerate(NWK):
        (cd / 'coassoc' / f'cover{c}.nwk').write_text(nwk + '\n')
        (cd / 'coassoc' / f'cover{c}.shards.tsv').write_text(f'c{c}.0000\t{len(TIPS)}\t0\n')
        (cd / 'trees' / f'c{c}.0000.trim.fasta').write_text('>x\n' + 'A' * 500 + '\n')
    return tmp_path


def test_coarse_group_spread_is_the_local_patristic_spread(tmp_path):
    rd = covers(tmp_path)
    h = {'tips': [{'id': t} for t in TIPS]}
    layer = {'unit_of_tip': [0, 0, 0, 0, 1, 1, 1], 'levels': [{'level': 0.125, 'groups': [{'units': [0, 1]}]}]}
    g = sparse.add_spread(rd, h, layer, procs=1, log=lambda *a: None)['levels'][0]['groups'][0]
    med = {}
    for a, b in itertools.combinations(TIPS, 2):
        d = [Tree(n, format=1).get_distance(a, b) for n in NWK]
        med[a, b] = (d[0] + d[1]) / 2                 # median of two trees
    assert abs(g['spread'] - max(med.values())) < 1e-4
    assert abs(g['spread_median'] - statistics.median(med.values())) < 1e-4
    assert g['spread_median'] < g['spread'] / 2       # most pairs are close; one far member sets the width


def test_a_layer_without_groups_is_left_as_it_is(tmp_path):
    layer = {'unit_of_tip': [0, 1], 'levels': [{'level': 0.125, 'groups': []}]}
    assert sparse.add_spread(tmp_path, {'tips': [{'id': 'a'}, {'id': 'b'}]}, layer) == layer
