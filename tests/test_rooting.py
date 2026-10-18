"""Midpoint rooting of shard trees (2026-10-17). A near-identical cluster (a1-a4) and a distant
relative (b1-b3), with the only scaffold tip (s) beside a1, as when a large component's scaffold is a
single one of its own centroids: rooted on the scaffold, the root falls inside the cluster and its
pairs join at nearly the whole tree's width; rooted at the midpoint they join at their own width."""
from seqhesion import coassoc, shardtrees

NWK = '((a1:0.002,s:0.001):0.003,a2:0.004,a3:0.003,a4:0.004,(b1:0.01,b2:0.012,b3:0.011):0.06);'
MEMBERS = ['a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'b3']


def joins(tmp_path, rooting):
    (tmp_path / 't.nwk').write_text(NWK)
    (tmp_path / 't.trim.fasta').write_text('>x\n' + 'A' * 500 + '\n')
    job = ('t', tmp_path / 't.nwk', tmp_path / 't.trim.fasta', MEMBERS, ['s'], rooting)
    _, line, n, _ = shardtrees.prepare(job)
    assert n == len(MEMBERS) and 's' not in line
    index = {t: i for i, t in enumerate(MEMBERS)}
    coassoc._init(index, 'identical')
    coassoc._G['floor'] = 1.0 / 500
    k, v = coassoc.joins(line)
    return {(MEMBERS[x // len(MEMBERS)], MEMBERS[x % len(MEMBERS)]): float(y) for x, y in zip(k, v)}


def test_scaffold_rooting_inflates_the_cluster(tmp_path):
    j = joins(tmp_path, 'scaffold')
    assert j[('a2', 'a3')] > 0.05                      # the whole tree's width


def test_midpoint_rooting_keeps_the_cluster_tight(tmp_path):
    j = joins(tmp_path, 'midpoint')
    inside = [v for (a, b), v in j.items() if a[0] == b[0] == 'a']
    across = [v for (a, b), v in j.items() if a[0] != b[0]]
    assert max(inside) < 0.02 < min(across)


def test_midpoint_is_the_default():
    assert shardtrees.ROOTING == 'midpoint'


def test_two_member_tree_is_left_as_is(tmp_path):
    (tmp_path / 't.nwk').write_text('(a1:0.01,a2:0.02,s:0.5);')
    (tmp_path / 't.trim.fasta').write_text('>x\n' + 'A' * 500 + '\n')
    _, line, n, _ = shardtrees.prepare(('t', tmp_path / 't.nwk', tmp_path / 't.trim.fasta', ['a1', 'a2'], ['s']))
    from ete3 import Tree
    assert n == 2 and sorted(Tree(line, format=5).get_leaf_names()) == ['a1', 'a2']
