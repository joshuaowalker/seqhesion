"""Coarse levels above the fine hierarchy: a layer over its units, from sparse shards (Josh, 2026-10-05).

Local shards (a seed and its nearest neighbours) almost never hold two distant tips, so above
~0.1 the fine hierarchy rests on the few near pairs they do hold, and average linkage over those
chains (lab: 5 groups at 0.2 in the 8,806-tip Cortinariaceae component). Here:

  units        the fine hierarchy's groups at UNIT_LEVEL (0.1), a lone tip being a unit of one
  sparse shard one random member of each of SIZE units drawn at random (stratified: every unit
               pair is equally likely to meet), plus a fixed outgroup from outside the component as
               rooting candidates, pruned like the scaffold; enough shards that each unit pair is
               expected in PAIR_VOTES of them. Trees: the production pipeline (cover.build_trees,
               shardtrees.export)
  layer        every shard's join level for each pair of its tips in different units is one vote
               on that unit pair; the median over its votes is the unit pair's level; average
               linkage over observed unit pairs gives a hierarchy of units, cut at COARSE_LEVELS.
               Unit pairs, not tip pairs, are weighed equally: feeding the sparse votes into the
               tip-level hierarchy instead let the near pairs local shards hold dominate, and its
               coarse levels chained (lab, experiments/seqhesion/merged_hierarchy.py)
  evidence     per coarse group of 2+ units: cohesion, pull, margin over unit-pair votes, as the
               fine levels over tip-pair votes (`pooled`, `group_pull`); held = share of its unit
               pairs observed. A fine group that is still a single unit at a coarse level is the
               same group (node) continuing; its internal evidence stays its fine-level one

Small components (under one shard): every local shard already holds all their own tips, so every
unit pair is already observed without bias; their layer is built from the local trees and they
get no sparse shards. Lab results (2026-10-05, main Cortinariaceae component, 102 trees): two
independent draws agree to 0.3-0.4 (ARI 0.87 at 0.3); random rather than stratified draws
observed a third of the unit pairs. Labels play no part.
"""
import collections
import hashlib
import json
import math
import random
import subprocess
from pathlib import Path

import numpy as np

from . import cover, shardtrees
from .coassoc import consensus_levels, hierarchy as linkage, pooled, votes_at
from .fasta import read_fasta, write_fasta
from .hierarchy import cut, group_pull, rnd

UNIT_LEVEL = 0.1
COARSE_LEVELS = (0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
SIZE = 150
PAIR_VOTES = 24
N_OUTGROUP = 8
METHOD = ('coarse layer: units = fine groups (and lone tips) at 0.1; stratified sparse shards (one random member '
          'of 150 random units, outside outgroup, ~24 votes per unit pair; small components: their local trees); '
          'unit-pair median join level, average linkage over observed unit pairs; v1 (2026-10-05)')


def units(h, level=UNIT_LEVEL):
    """(unit per tip index, {unit: [tip indices]}) at `level` of a hierarchy: its groups there,
    then each tip in none as a unit of one."""
    lv = next(x for x in h['levels'] if abs(x['level'] - level) < 1e-9)
    g = np.array(lv['group_of'])
    G = int(g.max()) + 1 if len(g) and g.max() >= 0 else 0
    unit = np.where(g >= 0, g, G + np.cumsum(g < 0) - 1)
    members = collections.defaultdict(list)
    for i, u in enumerate(unit.tolist()):
        members[u].append(i)
    return unit, members


def _seed(anchor):
    return int(hashlib.sha256(f'sparse:{anchor}'.encode()).hexdigest()[:12], 16)


def outgroups(regions, intake_dir, threads=8, log=print):
    """{region: [N_OUTGROUP tip ids]}: per region, the corpus sequences OUTSIDE its component hit
    most often (vsearch >= 0.5, best 400 per query) by its 85% centroids. One search for all."""
    d = Path(intake_dir)
    q, hits = d / 'sparse_outgroup_queries.fasta', d / 'sparse_outgroup_hits.tsv'
    own, queries = {}, {}
    for rd in regions:
        rd = Path(rd)
        own[rd.name] = set(read_fasta(rd / 'comp.fasta'))
        for h, s in read_fasta(rd / 'comp_cent_0.85.fasta').items():
            queries[f'{rd.name}|{h.split(";")[0]}'] = s
    write_fasta(q, queries)
    subprocess.run(['vsearch', '--usearch_global', str(q), '--db', str(d / 'full.derep.fasta'), '--id', '0.5',
                    '--maxaccepts', '400', '--maxrejects', '256', '--threads', str(threads), '--userout', f'{hits}.tmp',
                    '--userfields', 'query+target+id', '--quiet'], check=True, stderr=subprocess.DEVNULL)
    count = collections.defaultdict(collections.Counter)
    for line in open(f'{hits}.tmp'):
        qq, t, _ = line.rstrip('\n').split('\t')
        a = qq.split('|')[0]
        if t not in own[a]:
            count[a][t] += 1
    q.unlink()
    Path(f'{hits}.tmp').unlink()
    out = {}
    for rd in regions:
        a = Path(rd).name
        out[a] = [t for t, _ in sorted(count[a].items(), key=lambda x: (-x[1], x[0]))[:N_OUTGROUP]]
        if not out[a]:
            log(f'{a}: no sequence outside the component within 50% identity: rooted without an outgroup')
    return out


def plan(region_dir, h, outgroup, outgroup_seqs, size=SIZE, pair_votes=PAIR_VOTES):
    """Write REGION/sparse/{manifest.json, outgroup.fasta, units.json} for a large component.
    Returns the number of shards (0: fewer than 3 units, nothing written)."""
    rd = Path(region_dir)
    tips = [t['id'] for t in h['tips']]
    unit, members = units(h)
    U = len(members)
    if U < 3:
        return 0
    s = min(size, U)
    K = math.ceil(pair_votes / ((s / U) * ((s - 1) / (U - 1))))
    rng = random.Random(_seed(rd.name))
    shards = {}
    for k in range(K):
        pick = [rng.choice(members[u]) for u in rng.sample(range(U), s)]
        shards[f'sp{k:04d}'] = {'cover': 0, 'members': sorted(tips[i] for i in pick)}
    sd = rd / 'sparse'
    sd.mkdir(exist_ok=True)
    write_fasta(sd / 'outgroup.fasta', {t: outgroup_seqs[t] for t in outgroup}, outgroup)
    json.dump({'unit_level': UNIT_LEVEL, 'unit_of_tip': unit.tolist(), 'tips': tips}, open(sd / 'units.json', 'w'))
    json.dump({'scaffold': list(outgroup), 'shards': shards,
               'sparse': {'units': U, 'size': s, 'pair_votes': pair_votes, 'shards': K, 'seed': f'sparse:{rd.name}'}},
              open(sd / 'manifest.json', 'w'))
    return K


def build(region_dir, procs=8, log=print):
    """Trees for a planned region's sparse shards (REGION/sparse/trees, coassoc/), then its layer."""
    rd = Path(region_dir)
    sd = rd / 'sparse'
    man = json.load(open(sd / 'manifest.json'))
    seqs = read_fasta(rd / 'comp.fasta')
    seqs.update(read_fasta(sd / 'outgroup.fasta'))
    cover.build_trees(man, seqs, sd / 'trees', procs, log=log)
    shardtrees.export(sd, 0, procs, log=log)
    u = json.load(open(sd / 'units.json'))
    _, lines, widths = shardtrees.read(sd, 0)
    layer = build_layer(lines, widths, u['tips'], np.array(u['unit_of_tip']), log=log, procs=procs)
    layer['source'] = f'{len(lines)} sparse shard trees'
    _write(sd / 'layer.json', layer)
    return layer


def build_local(region_dir, h, covers='covers', procs=8, log=print):
    """The layer of a small component from its local trees (both covers, cross-cover twins once)."""
    rd = Path(region_dir)
    cd = rd / covers
    skip = cover.twins(json.load(open(cd / 'manifest.json')))
    lines, widths = [], []
    for c in (0, 1):
        for name, line, w in zip(*shardtrees.read(cd, c)):
            if name not in skip:
                lines.append(line)
                widths.append(w)
    tips = [t['id'] for t in h['tips']]
    unit, _ = units(h)
    layer = build_layer(lines, widths, tips, unit, log=log, procs=procs)
    layer['source'] = f'{len(lines)} local shard trees (small component)'
    (rd / 'sparse').mkdir(exist_ok=True)
    _write(rd / 'sparse' / 'layer.json', layer)
    return layer


def _write(path, obj):
    json.dump(obj, open(f'{path}.tmp', 'w'), separators=(',', ':'))
    Path(f'{path}.tmp').replace(path)


def build_layer(lines, widths, tips, unit, join='identical', procs=8, log=print):
    """The layer over units from prepared shard trees: unit-pair medians of every vote on a tip
    pair spanning two units, average linkage over observed unit pairs, cut at COARSE_LEVELS."""
    n, U = len(tips), int(unit.max()) + 1
    index = {t: i for i, t in enumerate(tips)}
    KK, VV = consensus_levels(lines, widths, index, procs, join)[3]
    ui, uj = unit[KK // n], unit[KK % n]
    keep = ui != uj
    UK = np.minimum(ui, uj)[keep].astype(np.int64) * U + np.maximum(ui, uj)[keep]
    UV = VV[keep]
    order = np.lexsort((UV, UK))
    UK, UV = UK[order], UV[order]
    starts = np.flatnonzero(np.r_[True, UK[1:] != UK[:-1]])
    counts = np.diff(np.r_[starts, len(UK)])
    ukeys = UK[starts]
    umed = (UV[starts + (counts - 1) // 2] + UV[starts + counts // 2]) / 2
    Z, observed, possible = linkage(ukeys, umed, U)
    levels = []
    for L in COARSE_LEVELS:
        cg = cut(Z, L)
        CG = int(cg.max()) + 1
        st = pooled(*votes_at(UK, UV, L), U, cg)
        pull = group_pull(st, cg, CG)
        groups = []
        for k in range(CG):
            us = np.flatnonzero(cg == k)
            coh = rnd(st['cohesion'][k])
            pairs = len(us) * (len(us) - 1) // 2
            groups.append({'units': us.tolist(), 'cohesion': coh, 'votes': int(st['c_votes'][k]), 'pull': rnd(pull[k]),
                           'margin': rnd(coh - pull[k]) if coh is not None else None,
                           'held': round(float(st['c_pairs'][k] / pairs), 3)})
        # per unit: its outside pull at this level ([share, coarse group or -1 for a lone unit, votes, the unit])
        unit_pull = [None if g == -2 else [rnd(v), int(g), int(pv), int(pt)]
                     for v, g, pv, pt in zip(st['pull'], st['p_group'], st['p_votes'], st['p_tip'])]
        levels.append({'level': L, 'unit_group': cg.tolist(), 'groups': groups, 'unit_pull': unit_pull,
                       'unit_membership': [rnd(x) for x in st['membership']], 'unit_m_votes': st['m_votes'].tolist()})
    log(f'layer: {U} units, {len(ukeys)} unit pairs observed ({len(ukeys) / max(1, U * (U - 1) // 2):.1%}), '
        f'median {int(np.median(counts)) if len(counts) else 0} votes; groups of 2+ units per level '
        + ', '.join(f'{lv["level"]}: {len(lv["groups"])}' for lv in levels))
    return {'method': METHOD, 'unit_level': UNIT_LEVEL, 'units': U, 'unit_of_tip': unit.tolist(),
            'unit_pairs_observed': len(ukeys), 'unit_pairs': U * (U - 1) // 2,
            'linkage': [[int(r[0]), int(r[1]), None if not np.isfinite(r[2]) else float(r[2]), int(r[3])] for r in Z],
            'observed': [int(x) for x in observed], 'possible': [int(x) for x in possible],
            'levels': levels}


def add_layer(h, layer):
    """The hierarchy `h` (hierarchy.build or its JSON) with the layer's coarse levels appended, in
    the same shape, so a release treats them like the fine ones. Per coarse level: groups of 2+
    units (layer evidence), plus every unit of 2+ tips still alone (the same node as its fine group,
    continuing); per tip, its unit's membership and pull. A node continuing from the fine levels
    keeps its fine measurement (kmax) and gains kmax_top, its coarsest level; nodes born in the
    layer carry no tree-level evidence (stem, spread, nearest, verdicts: None)."""
    tips = [t['id'] for t in h['tips']]
    assert layer['unit_level'] == UNIT_LEVEL and len(layer['unit_of_tip']) == len(tips)
    unit = np.array(layer['unit_of_tip'])
    U = layer['units']
    members = collections.defaultdict(list)
    for i, u in enumerate(unit.tolist()):
        members[u].append(i)
    node_of = {tuple(nd_tips): k for k, nd_tips in _node_tips(h).items()}
    fine_k = len(h['levels']) - 1
    for lv in h['levels']:
        lv.setdefault('layer', False)
    for L in layer['levels']:
        cg = np.array(L['unit_group'])
        CG = int(cg.max()) + 1 if len(cg) and cg.max() >= 0 else 0
        # coarse groups: layer groups first (2+ units), then lone units of 2+ tips, largest first
        sets = [(sorted(i for u in g['units'] for i in members[u]), g) for g in L['groups']]
        sets += [(sorted(members[u]), None) for u in range(U) if cg[u] < 0 and len(members[u]) >= 2]
        order = sorted(range(len(sets)), key=lambda k: (-len(sets[k][0]), sets[k][0][0]))
        group_of = np.full(len(tips), -1)
        unit_group = {}
        groups = []
        k_level = len(h['levels'])
        for gid, k in enumerate(order):
            tt, lg = sets[k]
            group_of[tt] = gid
            for u in {int(unit[i]) for i in tt}:
                unit_group[u] = gid
            key = tuple(tt)
            if key in node_of:
                nd = node_of[key]
                h['nodes'][nd]['kmax_top'] = k_level
            else:
                nd = len(h['nodes'])
                node_of[key] = nd
                h['nodes'].append({'kmin': k_level, 'kmax': k_level, 'layer': True, 'shards': [], 'verdicts': None,
                                   'stem': None, 'stem_n': 0, 'spread': None, 'spread_median': None, 'nearest': None,
                                   'second': None, 'stray': None, 'strays_closer': None, 'gap': None})
            groups.append({'id': gid, 'tips': tt, 'node': nd,
                           'cohesion': lg['cohesion'] if lg else None, 'votes': lg['votes'] if lg else 0,
                           'pull': lg['pull'] if lg else None, 'margin': lg['margin'] if lg else None,
                           'held': lg['held'] if lg else None})
        # a lone unit's group gets its own pull at this level
        for g in groups:
            if g['pull'] is None:
                p = L['unit_pull'][int(unit[g['tips'][0]])]
                g['pull'] = p[0] if p else 0.0

        def target(p):
            if p is None:
                return None
            share, tg, votes, pt = p
            if tg >= 0:
                tgt = unit_group[int(np.flatnonzero(cg == tg)[0])]
                return [share, tgt, votes, 0]
            if tg == -1:                                  # a lone unit: its group, or its single tip
                if pt in unit_group:
                    return [share, unit_group[pt], votes, 0]
                return [share, -1, votes, 0, int(members[pt][0])]
            return [share, tg, votes, 0]
        h['levels'].append({
            'level': L['level'], 'layer': True,
            'summary': {'groups': len(groups), 'tips_in_groups': int(np.sum(group_of >= 0)),
                        'largest': max((len(g['tips']) for g in groups), default=0), 'held_p10': None},
            'groups': groups, 'group_of': group_of.tolist(),
            'membership': [L['unit_membership'][u] for u in unit.tolist()],
            'm_votes': [L['unit_m_votes'][u] for u in unit.tolist()],
            'm_partners': [None] * len(tips),
            'pull': [target(L['unit_pull'][u]) for u in unit.tolist()],
        })
    # nesting over all levels again
    K = len(h['levels'])
    for k in range(fine_k, K):
        for grp in h['levels'][k]['groups']:
            t0 = grp['tips'][0]
            grp['parent'] = h['levels'][k + 1]['group_of'][t0] if k + 1 < K else None
            finer = h['levels'][k - 1]['group_of']
            grp['children'] = sorted({finer[t] for t in grp['tips'] if finer[t] >= 0})
    h['layer'] = {k: v for k, v in layer.items() if k in ('method', 'unit_level', 'units', 'unit_pairs_observed', 'unit_pairs', 'source')}
    h['linkage'] = combined_linkage(h['linkage'], layer, len(tips))
    return h


def combined_linkage(fine, layer, n):
    """The dendrogram every level is cut from: the fine linkage inside each unit (its subtree),
    the layer's linkage over units above them. A layer merge never sits below its children, nor at
    or below UNIT_LEVEL (it joins units the fine hierarchy keeps apart there): such heights are
    raised to just above it. Merge rows' observed / possible count unit pairs above the units."""
    unit = layer['unit_of_tip']
    Z, obs, pos = fine['Z'], fine['observed'], fine['possible']
    size_unit = collections.Counter(unit)
    rows, robs, rpos = [], [], []
    new_id, height, unit_of_node = {}, {}, {}
    for i in range(n):
        new_id[i], height[i], unit_of_node[i] = i, 0.0, unit[i]
    root_of_unit = {unit[i]: i for i in range(n) if size_unit[unit[i]] == 1}
    count = {i: 1 for i in range(n)}
    for r, (a, b, ht, sz) in enumerate(Z):
        a, b = int(a), int(b)
        ua, ub = unit_of_node.get(a), unit_of_node.get(b)
        if ua is None or ua != ub:
            unit_of_node[n + r] = None                   # a merge across units: replaced by the layer
            continue
        unit_of_node[n + r] = ua
        k = n + len(rows)
        rows.append([new_id[a], new_id[b], ht, count[a] + count[b]])
        robs.append(obs[r])
        rpos.append(pos[r])
        new_id[n + r], height[k], count[n + r] = k, ht or 0.0, count[a] + count[b]
        if count[n + r] == size_unit[ua]:
            root_of_unit[ua] = k
    floor = UNIT_LEVEL + 1e-6
    node = dict(root_of_unit)                           # layer node id -> combined node id
    size = {root_of_unit[u]: size_unit[u] for u in root_of_unit}
    hgt = {root_of_unit[u]: height.get(root_of_unit[u], 0.0) for u in root_of_unit}
    U = layer['units']
    for r, (a, b, ht, _) in enumerate(layer['linkage']):
        na, nb = node[int(a)], node[int(b)]
        k = n + len(rows)
        h_ = None if ht is None else max(ht, hgt[na], hgt[nb], floor)
        rows.append([na, nb, h_, size[na] + size[nb]])
        robs.append(layer['observed'][r])
        rpos.append(layer['possible'][r])
        node[U + r], size[k], hgt[k] = k, size[na] + size[nb], (h_ if h_ is not None else max(hgt[na], hgt[nb]))
    assert len(rows) == n - 1 or not layer['linkage'], 'combined linkage is not a full tree'
    return {'Z': rows, 'observed': robs, 'possible': rpos}


def _node_tips(h):
    out = {}
    for lv in h['levels']:
        for g in lv['groups']:
            out[g['node']] = g['tips']
    return out


def empty_layer(h):
    """The layer of a component with fewer than two units: every unit alone at every coarse level."""
    unit, members = units(h)
    U = len(members)
    lone = [-1] * U
    levels = [{'level': L, 'unit_group': lone, 'groups': [], 'unit_pull': [None] * U,
               'unit_membership': [None] * U, 'unit_m_votes': [0] * U} for L in COARSE_LEVELS]
    return {'method': METHOD, 'unit_level': UNIT_LEVEL, 'units': U, 'unit_of_tip': unit.tolist(),
            'unit_pairs_observed': 0, 'unit_pairs': 0, 'linkage': [], 'observed': [], 'possible': [], 'levels': levels,
            'source': 'fewer than two units'}
