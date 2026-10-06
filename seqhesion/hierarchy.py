"""The co-association hierarchy of one region, with its evidence.

Method: for every pair of tips, each shard holding both records the diameter of the smallest
clade holding them (inside a polytomy, tips closer than one expected change -- which the data
cannot tell apart -- join at their own span: `coassoc.joins`, 'identical'); the median over
shards is the pair's join level; average linkage over observed pairs only builds one hierarchy,
and parts no observed pair connects stay a forest. It is cut at a series of levels.

The shards are POOLED over both covers (Josh, 2026-10-04): every distinct shard tree counts once.
The two covers are drawn independently, so dense neighbourhoods often give both the same shard
(55-61% of cover 1 over the Agaricales build); such a twin (`cover.twins`) is the same sample and
the same tree, counted once. Until v6 the hierarchy used cover 0 alone and cover 1 'replicated'
it; that one fixed split-half agreement was dropped: on 334 components it did not predict which
groups survive resampling the shards, while the margin below did (lab experiments, 2026-10-04).

All shares are POOLED: every shard's verdict on a pair is one vote (`coassoc.pooled`).
Per group and level:
  cohesion      share of the votes on the group's pairs that join the pair at or below the
                level, with the votes behind it
  pull          the members' outside pulls, POOLED like cohesion: over every member, the votes
                joining it to its strongest outside target (tip pull below) at or below the level,
                over the votes on those pairs (`group_pull`). Not the maximum over members: a
                member's pull often rests on a handful of votes (3 of 3 = 1.0), so a maximum
                saturated in large groups (mm-to-ref found it in release 20261004.03f)
  margin        cohesion - pull: how much more the shards join the group's own pairs than they
                join its members to anything outside. The confidence measure: it predicts which
                groups recur when the shards are resampled (AUC 0.94 for recurring in under half
                of 20 bootstrap resamples, 0.87 for under 95%; 0.83 in groups of 101+ tips, where
                the maximum managed 0.69) better than cohesion (0.80)
  held          share of the group's pairs that some shard holds
Per tip and level:
  membership    the same over its pairs with the rest of its group, with votes and partners
  pull          the same against the outside group (or ungrouped tip) it scores highest with
Per distinct group (a node: the same tips over a run of levels), from the pooled shard trees:
  shards        every shard holding >= 2 of its tips: tips held, outsiders held, verdict (clade /
                unresolved / conflict / whole), stem when a clade (expected changes), and the
                pairs it joins at the node's coarsest level
  spread        widest (and median) distance between its tips; distances are the median over
                shards of each pair's tip-to-tip distance in the shard trees
  nearest       the nearest outside GROUP (2+ tips) at the node's coarsest level, via its closest
                tip (so a subgroup is never a relative), and the second; gap = nearest / spread.
                Ungrouped tips are reported apart: the nearest one, and how many sit closer.

Stem nodes (adopted 2026-10-06, with mm-to-ref): a node of the linkage at height <= the top fine
level that is not a group at any level (it lives between two levels) but whose linkage stem (its
parent's height - its own) is >= STEM_MIN. A clade can live wholly between two levels (T. versicolor:
0.021-0.026); such a node is reported with the same cohesion / pull / margin / held, measured at the
largest GRID point below its parent's height (always inside its life, as STEM_MIN >= 2 GRID).
Under bootstrap resampling of the shard trees, 75% of stem-node ids carried against 90% of level
groups'; a longer stem bought no stability, so the threshold is the smallest the grid allows
(lab experiment, 2026-10-06). No shard verdicts, spread, nearest or identity are computed for them.

Labels play no part: tips are ids.
"""
import collections
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import cover, provenance, shardtrees
from .coassoc import VERDICTS, consensus_levels, hierarchy, pooled, separation, shard_evidence, votes_at

METHOD = ('co-association hierarchy: per pair, median over shards holding both of the diameter of '
          'the smallest clade holding them, except that inside a polytomy tips closer than one expected '
          'change join at their own span (join "identical"); average linkage over observed pairs only; shares pooled '
          'over shard votes; per node: shard evidence, spread and nearest relatives from median '
          'patristic distances, nearest relative = nearest outside group; stand-alone (no stitch); shards pooled over '
          'both covers, cross-cover twins once; per group: cohesion, pull (members\' outside pulls pooled over votes), '
          'margin = cohesion - pull (no replication); v7 (2026-10-04)')

LEVELS = (0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1)
STEM_MIN = 0.001                 # a stem node's linkage stem is at least this
GRID = 0.0005                    # stem nodes are measured at a multiple of this

_W = {}


def _init(index, width, groups_of_tip, group_mask, group_level, levels, join):
    _W.update(index=index, width=width, groups_of_tip=groups_of_tip, group_mask=group_mask,
              group_level=group_level, levels=levels, join=join)


def _evidence(job):
    from ete3 import Tree
    k, line = job
    t = Tree(line, format=5)
    return k, shard_evidence(t, _W['index'], _W['width'][k], _W['groups_of_tip'], _W['group_mask'],
                             _W['group_level'], _W['levels'], _W['join'])


def rnd(x, d=3):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), d)


def group_pull(s, group_of, G):
    """Per group, its members' outside pulls pooled over votes (`coassoc.pooled` output `s`): the
    votes joining each member to its strongest outside target, over the votes on those pairs. A
    member with nothing outside joining counts its outside votes with no hits; one with no outside
    votes counts nothing. A group with no outside votes at all: 0."""
    v = np.where(s['p_group'] == -2, 0, s['p_votes']).astype(float)
    h = np.rint(np.nan_to_num(s['pull']) * v)
    inside = group_of >= 0
    H = np.bincount(group_of[inside], h[inside], G)
    V = np.bincount(group_of[inside], v[inside], G)
    return np.where(V > 0, H / np.maximum(V, 1), 0.0)


def cut(Z, level):
    """Groups of 2+ tips at a level, numbered largest first; -1 for a tip in none."""
    from scipy.cluster.hierarchy import fcluster
    p = fcluster(Z, t=level, criterion='distance')
    size = collections.Counter(p.tolist())
    order = sorted((g for g, s in size.items() if s >= 2), key=lambda g: (-size[g], g))
    gid = {g: k for k, g in enumerate(order)}
    return np.array([gid.get(g, -1) for g in p.tolist()])


def linkage_nodes(Z, n, top):
    """[(row, tip indices, height, parent height)] for every row of linkage Z at height <= top
    (rows of height None / inf join forest parts; a missing parent is inf)."""
    h = [np.inf if r[2] is None or not np.isfinite(r[2]) else float(r[2]) for r in Z]
    parent_h = [np.inf] * len(Z)
    for r, row in enumerate(Z):
        for c in (int(row[0]), int(row[1])):
            if c >= n:
                parent_h[c - n] = h[r]
    members, out = {}, []
    for r, row in enumerate(Z):
        a, b = int(row[0]), int(row[1])
        s = (members.pop(a) if a >= n else [a]) + (members.pop(b) if b >= n else [b])
        members[n + r] = s
        if h[r] <= top:
            out.append((r, sorted(s), h[r], parent_h[r]))
    return out


def stem_nodes(Z, n, kv, levels, stem_min=STEM_MIN):
    """The stem nodes of linkage Z (see the module docstring), each measured like a group: a list of
    {tips, height, top, measured_at, cohesion, votes, pull, margin, held}, tips as indices. kv: the
    per-shard join levels (consensus_levels' fourth value)."""
    assert stem_min >= 2 * GRID
    top_level = max(levels)
    at = collections.defaultdict(list)
    for _, s, ht, ph in linkage_nodes(Z, n, top_level):
        if ph - ht >= stem_min and not any(ht <= L < ph for L in levels):
            g = round(float(np.floor((ph - 1e-12) / GRID) * GRID), 6)
            assert ht <= g < ph, (ht, ph, g)
            at[g].append((s, ht, ph))
    Zf = np.array([[r[0], r[1], np.inf if r[2] is None else r[2], r[3]] for r in Z], float)
    out = []
    for g, todo in sorted(at.items()):
        group_of = cut(Zf, g)
        G = int(group_of.max()) + 1
        s0 = pooled(*votes_at(*kv, g), n, group_of)
        pull = group_pull(s0, group_of, G)
        for s, ht, ph in todo:
            k = int(group_of[s[0]])
            assert k >= 0 and np.flatnonzero(group_of == k).tolist() == s, ('stem node is not a group at', g)
            coh = rnd(s0['cohesion'][k])
            out.append({'tips': s, 'height': round(ht, 6), 'top': round(ph, 6), 'measured_at': g, 'cohesion': coh,
                        'votes': int(s0['c_votes'][k]), 'pull': rnd(pull[k]),
                        'margin': rnd(coh - pull[k]) if coh is not None else None,
                        'held': round(float(s0['c_pairs'][k] / (len(s) * (len(s) - 1) // 2)), 3)})
    return out


def pooled_trees(cd):
    """(manifest, cover-1 twins, shard names, tree lines, trimmed widths): every distinct prepared
    shard tree of both covers, a cover-1 twin of a cover-0 shard counted once."""
    man = json.load(open(cd / 'manifest.json'))
    skip = cover.twins(man)
    shard_names, lines, width = [], [], []
    for c in (0, 1):
        for name, line, w in zip(*shardtrees.read(cd, c)):
            if name not in skip:
                shard_names.append(name)
                lines.append(line)
                width.append(w)
    return man, skip, shard_names, lines, width


def add_stem_nodes(h, covers_dir, procs=8, log=print):
    """Stem nodes for a hierarchy built before they existed: the pooled join levels are recomputed
    from the stored shard trees (no tree is built) and must reproduce its linkage."""
    tips = [t['id'] for t in h['tips']]
    n, index = len(tips), {t: i for i, t in enumerate(tips)}
    _, _, _, lines, width = pooled_trees(Path(covers_dir))
    keys, med, _, kv = consensus_levels(lines, width, index, procs, h['join'])
    Z = hierarchy(keys, med, n)[0]
    same = lambda A: {tuple(s) for _, s, *_ in linkage_nodes(A, n, np.inf)}  # noqa: E731
    assert same([[r[0], r[1], None if not np.isfinite(r[2]) else r[2], r[3]] for r in Z]) == same(h['linkage']['Z']), \
        'the stored shard trees do not reproduce this hierarchy'
    h['stem_nodes'] = stem_nodes(h['linkage']['Z'], n, kv, [lv['level'] for lv in h['levels']])
    h['stem_min'] = STEM_MIN
    log(f'{len(h["stem_nodes"])} stem nodes')
    return h


def build(covers_dir, tips, levels=LEVELS, join='identical', procs=8, log=print):
    """The hierarchy over `tips` (ids, in index order) from a covers directory whose shard trees
    have been prepared (`shardtrees.export`, both covers), pooled: every distinct shard tree of both
    covers, a cover-1 twin of a cover-0 shard (`cover.twins`) counted once. Returns a dict; tips are
    referred to by index into `tips`, shards by index into `shards`."""
    from scipy.cluster.hierarchy import fcluster
    cd = Path(covers_dir)
    levels = sorted(float(x) for x in levels)
    n = len(tips)
    index = {t: i for i, t in enumerate(tips)}

    man, skip, shard_names, lines, width = pooled_trees(cd)
    keys, med, _, kv = consensus_levels(lines, width, index, procs, join)
    Z, observed, possible = hierarchy(keys, med, n)
    log(f'{len(lines)} distinct shard trees (both covers; {len(skip)} cover-1 twins counted once): '
        f'{len(keys):,} co-sampled pairs')
    dkeys, dist, _, _ = consensus_levels(lines, width, index, procs, 'patristic')
    # the forest: parts of the hierarchy that no observed pair connects (0 = the largest)
    top = fcluster(Z, t=np.finfo(float).max, criterion='distance')
    part_size = collections.Counter(top.tolist())
    part_rank = {c: r for r, (c, _) in enumerate(sorted(part_size.items(), key=lambda x: (-x[1], x[0])))}
    part = [part_rank[c] for c in top.tolist()]
    seen = np.zeros(n, bool)
    seen[np.unique(np.r_[keys // n, keys % n])] = True

    out_levels = []
    for L in levels:
        group_of = cut(Z, L)
        G = int(group_of.max()) + 1
        s0 = pooled(*votes_at(*kv, L), n, group_of)
        pull = group_pull(s0, group_of, G)
        groups = []
        for k in range(G):
            members = np.flatnonzero(group_of == k)
            pairs = len(members) * (len(members) - 1) // 2
            coh = rnd(s0['cohesion'][k])
            groups.append({'id': k, 'tips': members.tolist(),
                           'cohesion': coh, 'votes': int(s0['c_votes'][k]),
                           'pull': rnd(pull[k]),
                           'margin': rnd(coh - pull[k]) if coh is not None else None,
                           'held': round(float(s0['c_pairs'][k] / pairs), 3)})
        out_levels.append({
            'level': L,
            'summary': {'groups': G, 'tips_in_groups': int(np.sum(group_of >= 0)),
                        'largest': int(max((len(g['tips']) for g in groups), default=0)),
                        'held_p10': round(float(np.quantile([g['held'] for g in groups], 0.1)), 3) if groups else None},
            'groups': groups,
            'group_of': group_of.tolist(),
            'membership': [rnd(x) for x in s0['membership']],
            'm_votes': s0['m_votes'].tolist(), 'm_partners': s0['m_partners'].tolist(),
            # [share, target, votes, partners]: target >= 0 a group, -1 an ungrouped tip (5th: the tip),
            # -3 outside evidence with nothing joining at this level; None = no outside evidence
            'pull': [None if g == -2 else [round(float(v), 3), int(g), int(pv), int(pp)] + ([int(pt)] if g == -1 else [])
                     for v, g, pv, pp, pt in zip(s0['pull'], s0['p_group'], s0['p_votes'], s0['p_partners'], s0['p_tip'])],
        })
        log(f'level {L}: {G} groups over {int(np.sum(group_of >= 0))} tips')

    # nesting: each group's parent at the next coarser level, children at the next finer one
    K = len(levels)
    for k in range(K):
        for grp in out_levels[k]['groups']:
            t0 = grp['tips'][0]
            grp['parent'] = out_levels[k + 1]['group_of'][t0] if k + 1 < K else None
            finer = out_levels[k - 1]['group_of'] if k > 0 else None
            grp['children'] = sorted({finer[t] for t in grp['tips'] if finer[t] >= 0}) if k > 0 else []

    # nodes: a group with the same tips over a run of levels is one node; its coarsest level (kmax)
    # is the one it is measured at
    node_of, nodes = {}, []
    for k in range(K):
        for grp in out_levels[k]['groups']:
            key = tuple(grp['tips'])
            if key not in node_of:
                node_of[key] = len(nodes)
                nodes.append({'tips': grp['tips'], 'kmin': k, 'kmax': k})
            nd = nodes[node_of[key]]
            nd['kmax'] = k
            grp['node'] = node_of[key]
    log(f'{len(nodes)} distinct groups (nodes)')

    # shard evidence per node, from the pooled shard trees
    tip_shards = [[] for _ in range(n)]                  # which shards hold each tip
    for k, s_name in enumerate(shard_names):
        for t in man['shards'][s_name]['members']:
            tip_shards[index[t]].append(k)
    groups_of_tip = [[] for _ in range(n)]
    for g, nd in enumerate(nodes):
        for t in nd['tips']:
            groups_of_tip[t].append(g)
    group_mask = [sum(1 << t for t in nd['tips']) for nd in nodes]
    group_level = [nd['kmax'] for nd in nodes]
    rows = [[] for _ in nodes]
    with Pool(procs, initializer=_init, initargs=(index, width, groups_of_tip, group_mask, group_level, levels, join)) as pool:
        for k, res in pool.imap_unordered(_evidence, list(enumerate(lines)), chunksize=8):
            for g, held, out, verdict, stem, hits, votes in res:
                rows[g].append([k, held, out, VERDICTS.index(verdict), rnd(stem, 2), hits, votes])
    log('shard evidence gathered')
    floor = 1.0 / float(np.median(width))           # about one expected change, in subs/site
    for g, nd in enumerate(nodes):
        r = sorted(rows[g], key=lambda x: (x[3], -x[1], x[0]))
        stems = [x[4] for x in r if x[3] == 0 and x[4] is not None]
        sep = separation(dkeys, dist, n, nd['tips'], np.array(out_levels[nd['kmax']]['group_of']))
        near = sep['nearest']
        nd.update({
            'shards': r,
            'verdicts': [sum(1 for x in r if x[3] == v) for v in range(len(VERDICTS))],
            'stem': rnd(float(np.median(stems)), 2) if stems else None, 'stem_n': len(stems),
            'spread': rnd(sep['spread'], 4), 'spread_median': rnd(sep['spread_median'], 4),
            'nearest': [rnd(near[0], 4), near[1], near[2]] if near else None,
            'second': [rnd(sep['second'][0], 4), sep['second'][1], sep['second'][2]] if sep['second'] else None,
            'stray': [rnd(sep['stray'][0], 4), sep['stray'][1]] if sep['stray'] else None,
            'strays_closer': sep['strays_closer'],
            'gap': rnd(near[0] / max(sep['spread'] or 0.0, floor), 2) if near else None,
        })
        del nd['tips']

    Zj = [[int(r[0]), int(r[1]), (None if not np.isfinite(r[2]) else float(r[2])), int(r[3])] for r in Z]
    stems = stem_nodes(Zj, n, kv, levels)
    log(f'{len(stems)} stem nodes')

    return {
        'method': METHOD, 'covers': 'pooled', 'twins': len(skip), 'join': join,
        'stem_nodes': stems, 'stem_min': STEM_MIN,
        'provenance': provenance.hierarchy_key(cd, METHOD, {'levels': levels, 'join': join, 'distance': 'patristic median'}),
        'verdicts': list(VERDICTS), 'shards': shard_names, 'distance_floor': round(floor, 5),
        'forest': sorted(part_size.values(), reverse=True),
        'tips': [{'id': t, 'seen': bool(seen[i]), 'part': part[i], 'shards': tip_shards[i]} for i, t in enumerate(tips)],
        'levels': out_levels,
        'nodes': nodes,
        'shard_trees': lines,
        # the dendrogram the levels are cut from: scipy linkage rows [a, b, height,
        # size], height None for rows joining forest parts (no observed pair between them), with the
        # observed and possible cross pairs behind each merge
        'linkage': {'Z': [[int(r[0]), int(r[1]), (None if not np.isfinite(r[2]) else float(r[2])), int(r[3])]
                          for r in Z],
                    'observed': [int(x) for x in observed],
                    'possible': [int(x) for x in possible]},
    }
