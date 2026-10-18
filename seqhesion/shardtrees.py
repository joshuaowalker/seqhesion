"""Each shard tree as co-association reads it.

  * scaffold and context tips PRUNED: they sit in every shard and would link to everything;
  * local shards ROOTED at the midpoint of the members' own tree (ROOTING, adopted 2026-10-17). Before, the tree
    was rooted on the scaffold tip farthest from the first 40 members, then pruned ('scaffold').
    For large components the scaffold is the component's own 85% centroids, often one or two:
    the root then fell inside a near-identical cluster, which became a grade with its relatives
    nested inside, and its pairs joined at the whole tree's width (lab, Hygrocybe 'miniata-IN01':
    0.124 vs 0.030 midpoint-rooted). Outgroups from outside the component failed differently: a
    distant outgroup attached inside a dense cluster in every tree (lab, a 775-tip Trametes
    component). The midpoint fixed both and left components already rooted well unchanged. Sparse and
    corpus shards (one member per unit) keep 'scaffold' rooting on their outside outgroups;
  * internal branches shorter than 0.5 / alignment columns (about half an expected change)
    COLLAPSED into polytomies, branch length kept on the children.

One Newick per line in <covers>/coassoc/cover<k>.nwk, and a TSV naming each line's shard
(name, member tips, branches collapsed).
"""
import json
from multiprocessing import Pool
from pathlib import Path

from .tools import alignment_width


ROOTING = 'midpoint'          # or 'scaffold' (until 2026-10-17), kept for comparison


def prepare(job):
    """(name, prepared Newick, member tips, branches collapsed) for one shard tree. job: (name, tree path,
    trimmed alignment path, members, scaffold + context tips[, rooting])."""
    from ete3 import Tree
    name, tree_path, trim_path, members, outgroup = job[:5]
    rooting = job[5] if len(job) > 5 else ROOTING
    tree = Tree(str(tree_path), format=0)
    leaves = {l.name: l for l in tree.iter_leaves()}
    mem = [leaves[m] for m in members if m in leaves]
    if rooting == 'scaffold':
        scaf = [leaves[s] for s in outgroup if s in leaves]
        if scaf and mem:
            far = max(scaf, key=lambda s: sum(s.get_distance(m) for m in mem[:40]))
            tree.set_outgroup(far)
        tree.prune([m.name for m in mem], preserve_branch_length=True)
    elif rooting == 'midpoint':
        tree.prune([m.name for m in mem], preserve_branch_length=True)
        if len(mem) > 2:
            mid = tree.get_midpoint_outgroup()
            if mid is not None and mid is not tree:
                tree.set_outgroup(mid)
    else:
        raise ValueError(f'rooting {rooting!r}')
    min_len = 0.5 / alignment_width(trim_path)
    collapsed = 0
    for node in list(tree.traverse('postorder')):
        if not node.is_leaf() and not node.is_root() and node.dist < min_len:
            node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
            collapsed += 1
    return name, tree.write(format=5), len(mem), collapsed


def paths(covers_dir, cover):
    d = Path(covers_dir) / 'coassoc'
    return d / f'cover{cover}.nwk', d / f'cover{cover}.shards.tsv'


def export(covers_dir, cover, procs=8, log=print, rooting=None):
    """Prepare one cover's trees into <covers>/coassoc. rooting: ROOTING (midpoint) by default, for local
    shards; the sparse and corpus layers pass 'scaffold': their shards hold one member per unit, so no
    dense cluster for an outgroup to fall into, and their outside outgroups carry the deep order."""
    cd = Path(covers_dir)
    man = json.load(open(cd / 'manifest.json'))
    scaffold = man.get('scaffold', [])
    rooting = rooting or ROOTING
    jobs = [(n, cd / 'trees' / f'{n}.nwk', cd / 'trees' / f'{n}.trim.fasta', s['members'], scaffold + s.get('context', []),
             rooting) for n, s in sorted(man['shards'].items()) if s['cover'] == cover]
    nwk, idx = paths(cd, cover)
    nwk.parent.mkdir(exist_ok=True)
    tips, coll = 0, 0
    with Pool(procs) as pool, open(f'{nwk}.tmp', 'w') as f, open(f'{idx}.tmp', 'w') as g:
        for name, newick, n, c in pool.imap(prepare, jobs, chunksize=8):
            f.write(newick + '\n')
            g.write(f'{name}\t{n}\t{c}\n')
            tips += n
            coll += c
    Path(f'{nwk}.tmp').replace(nwk)
    Path(f'{idx}.tmp').replace(idx)
    log(f'wrote {nwk}: {len(jobs)} trees, {tips} member tips in all, {coll} short branches collapsed')


def read(covers_dir, cover):
    """(shard names, Newick lines, trimmed alignment widths) for one cover."""
    cd = Path(covers_dir)
    nwk, idx = paths(cd, cover)
    names = [l.split('\t')[0] for l in open(idx)]
    lines = [l for l in open(nwk) if l.strip()]
    assert len(names) == len(lines), f'{nwk}: {len(lines)} trees, {len(names)} names'
    widths = [alignment_width(cd / 'trees' / f'{s}.trim.fasta') for s in names]
    return names, lines, widths
