"""The co-association hierarchy of one region, with its evidence.

Method: for every pair of tips, each shard holding both records the diameter of the smallest
clade holding them (inside a polytomy, tips closer than one expected change -- which the data
cannot tell apart -- join at their own span: `coassoc.joins`, 'identical'); the median over
shards is the pair's join level; average linkage over observed pairs only builds one hierarchy
from the primary cover, and parts no observed pair connects stay a forest. It is cut at a series
of levels.

All shares are POOLED: every shard's verdict on a pair is one vote (`coassoc.pooled`).
Per group and level:
  cohesion      share of the votes on the group's pairs, in the primary cover, that join the
                pair at or below the level, with the votes behind it
  replication   the same over the other cover, which the hierarchy never saw
  held          share of the group's pairs that some primary-cover shard holds
Per tip and level:
  membership    the same over its pairs with the rest of its group, with votes and partners
  pull          the same against the outside group (or ungrouped tip) it scores highest with
Per distinct group (a node: the same tips over a run of levels), from the primary cover's trees:
  shards        every shard holding >= 2 of its tips: tips held, outsiders held, verdict (clade /
                unresolved / conflict / whole), stem when a clade (expected changes), and the
                pairs it joins at the node's coarsest level
  spread        widest (and median) distance between its tips; distances are the median over
                shards of each pair's tip-to-tip distance in the shard trees
  nearest       the nearest outside GROUP (2+ tips) at the node's coarsest level, via its closest
                tip (so a subgroup is never a relative), and the second; gap = nearest / spread.
                Ungrouped tips are reported apart: the nearest one, and how many sit closer.

Labels play no part: tips are ids.
"""
import collections
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import provenance, shardtrees
from .coassoc import VERDICTS, consensus_levels, hierarchy, pooled, separation, shard_evidence, votes_at
from .stats import adjusted_rand

METHOD = ('co-association hierarchy: per pair, median over shards holding both of the diameter of '
          'the smallest clade holding them, except that inside a polytomy tips closer than one expected '
          'change join at their own span (join "identical"); average linkage over observed pairs only; shares pooled '
          'over shard votes; per node: shard evidence, spread and nearest relatives from median '
          'patristic distances, nearest relative = nearest outside group; stand-alone (no stitch); v5 (2026-09-30)')

LEVELS = (0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1)

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


def cut(Z, level):
    """Groups of 2+ tips at a level, numbered largest first; -1 for a tip in none."""
    from scipy.cluster.hierarchy import fcluster
    p = fcluster(Z, t=level, criterion='distance')
    size = collections.Counter(p.tolist())
    order = sorted((g for g, s in size.items() if s >= 2), key=lambda g: (-size[g], g))
    gid = {g: k for k, g in enumerate(order)}
    return np.array([gid.get(g, -1) for g in p.tolist()])


def build(covers_dir, tips, levels=LEVELS, join='identical', primary=0, procs=8, log=print):
    """The hierarchy over `tips` (ids, in index order) from a covers directory whose shard trees
    have been prepared (`shardtrees.export`, both covers). Returns a dict; tips are referred to by
    index into `tips`, shards by index into `shards`."""
    from scipy.cluster.hierarchy import fcluster
    cd = Path(covers_dir)
    levels = sorted(float(x) for x in levels)
    n = len(tips)
    index = {t: i for i, t in enumerate(tips)}

    P, R = primary, 1 - primary          # primary and replicate cover
    read = {c: shardtrees.read(cd, c) for c in (0, 1)}
    covers = {}
    for c in (0, 1):
        _, lines, widths = read[c]
        keys, med, _, kv = consensus_levels(lines, widths, index, procs, join)
        Z, obs, pos = hierarchy(keys, med, n)
        covers[c] = {'keys': keys, 'kv': kv, 'Z': Z, 'observed': obs, 'possible': pos}
        log(f'cover {c}: {len(keys):,} co-sampled pairs')
    shard_names, lines, width = read[P]
    dkeys, dist, _, _ = consensus_levels(lines, width, index, procs, 'patristic')
    # the forest: parts of the primary hierarchy that no observed pair connects (0 = the largest)
    top = fcluster(covers[P]['Z'], t=np.finfo(float).max, criterion='distance')
    part_size = collections.Counter(top.tolist())
    part_rank = {c: r for r, (c, _) in enumerate(sorted(part_size.items(), key=lambda x: (-x[1], x[0])))}
    part = [part_rank[c] for c in top.tolist()]
    seen = [np.zeros(n, bool), np.zeros(n, bool)]
    for c in (0, 1):
        k = covers[c]['keys']
        seen[c][np.unique(np.r_[k // n, k % n])] = True
    both_seen = seen[0] & seen[1]

    out_levels = []
    for L in levels:
        lab = {c: cut(covers[c]['Z'], L) for c in (0, 1)}
        group_of = lab[P]
        G = int(group_of.max()) + 1
        s0 = pooled(*votes_at(*covers[P]['kv'], L), n, group_of)
        s1 = pooled(*votes_at(*covers[R]['kv'], L), n, group_of)     # replication: the other cover's votes
        groups = []
        for k in range(G):
            members = np.flatnonzero(group_of == k)
            pairs = len(members) * (len(members) - 1) // 2
            groups.append({'id': k, 'tips': members.tolist(),
                           'cohesion': rnd(s0['cohesion'][k]), 'votes': int(s0['c_votes'][k]),
                           'replication': rnd(s1['cohesion'][k]), 'votes1': int(s1['c_votes'][k]),
                           'seen1': round(float(s1['c_pairs'][k] / pairs), 3),
                           'held': round(float(s0['c_pairs'][k] / pairs), 3)})
        a0, a1 = np.where(group_of >= 0, group_of, -1 - np.arange(n)), np.where(lab[R] >= 0, lab[R], -1 - np.arange(n))
        out_levels.append({
            'level': L,
            'summary': {'groups': G, 'tips_in_groups': int(np.sum(group_of >= 0)),
                        'ari_covers': round(adjusted_rand(a0[both_seen], a1[both_seen]), 3),
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
        log(f'level {L}: {G} groups over {int(np.sum(group_of >= 0))} tips; ARI covers {out_levels[-1]["summary"]["ari_covers"]}')

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

    # shard evidence per node, from the primary cover's prepared trees
    man = json.load(open(cd / 'manifest.json'))
    assert all(man['shards'][s]['cover'] == P for s in shard_names)
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

    return {
        'method': METHOD, 'primary_cover': P, 'join': join,
        'provenance': provenance.hierarchy_key(cd, METHOD, {'levels': levels, 'join': join, 'distance': 'patristic median'}),
        'verdicts': list(VERDICTS), 'shards': shard_names, 'distance_floor': round(floor, 5),
        'forest': sorted(part_size.values(), reverse=True),
        'tips': [{'id': t, 'seen': bool(seen[P][i]), 'part': part[i], 'shards': tip_shards[i]} for i, t in enumerate(tips)],
        'levels': out_levels,
        'nodes': nodes,
        'shard_trees': lines,
        # the dendrogram the levels are cut from (primary cover): scipy linkage rows [a, b, height,
        # size], height None for rows joining forest parts (no observed pair between them), with the
        # observed and possible cross pairs behind each merge
        'linkage': {'Z': [[int(r[0]), int(r[1]), (None if not np.isfinite(r[2]) else float(r[2])), int(r[3])]
                          for r in covers[P]['Z']],
                    'observed': [int(x) for x in covers[P]['observed']],
                    'possible': [int(x) for x in covers[P]['possible']]},
    }
