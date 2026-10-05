"""Corpus levels: structure ACROSS components, from a cover over super-units (Josh, 2026-10-05).

Components are cut apart at 86% centroid identity, so nothing relates them. Here:

  super-units  every component's groups at UNIT_LEVEL (0.3, its coarse layer, seqhesion.sparse),
               a tip in none being a unit of one: ~4,800 over the corpus
  neighbours   up to REPS representatives per unit (its tips with most inputs); a unit's
               neighbours are the other units ranked by the best identity between representatives
               (vsearch; identity only chooses what is sampled together, never a grouping)
  cover        cover.plan over units as it plans tips: each shard a seed unit and its nearest
               SIZE - 1 units, COVERS covers of DEPTH, so units within reach meet in many shards
               and distant ones (ITS unalignable across orders) never do. A shard holds one random
               member of each of its units, plus one random member of each of the seed's next
               CONTEXT units as rooting candidates, pruned like the scaffold
  layer        as the component layer (sparse.build_layer): unit-pair median join levels, average
               linkage over observed unit pairs only (units nothing connects stay a forest), cut
               at LEVELS, with cohesion / pull / margin per group

Lab basis (experiments/seqhesion/cross_shards.py, 2026-10-05): over 122 components near the main
Cortinariaceae, two independent draws agreed at 0.3-0.75 (ARI 0.90-0.98) and not at 1.0 (0.66);
families (display only) came out as single nodes (best-node F1 0.96-0.97 vs 0.82 apart). So
LEVELS stop at 0.75. A navigation aid, not a phylogeny. Labels play no part.

Steps: plan (needs the corpus: run where the intake is), trees (chunked: --chunk I --of N, e.g. on
AWS), layer (prepares the trees and builds the layer).
"""
import collections
import json
import random
import subprocess
from pathlib import Path

import numpy as np

from . import cover, regions, shardtrees, sparse
from .fasta import read_fasta, write_fasta

UNIT_LEVEL = 0.3
LEVELS = (0.4, 0.5, 0.75)
REPS = 3
NBR_ID = 0.5
SIZE = 150
COVERS = 2
DEPTH = 6
CONTEXT = 20
METHOD = ('corpus layer: super-units = each component\'s groups (and lone tips) at 0.3; a cover over units by '
          'representative identity (2 covers, depth 6, 150 units per shard, one random member each, the seed\'s next 20 '
          'units as rooting context); unit-pair median join level, average linkage over observed unit pairs; v1 (2026-10-05)')


def super_units(region_dirs, hierarchy_name, log=print):
    """(tips, unit per tip, component per tip): each region's groups at UNIT_LEVEL after its coarse
    layer, then its tips in none, each a unit of one."""
    tips, unit, comp = [], [], []
    nxt = 0
    for rd in region_dirs:
        rd = Path(rd)
        h = sparse.add_layer(json.load(open(rd / hierarchy_name)), json.load(open(rd / 'sparse' / 'layer.json')))
        lv = next(x for x in h['levels'] if abs(x['level'] - UNIT_LEVEL) < 1e-9)
        g = np.array(lv['group_of'])
        G = int(g.max()) + 1 if (g >= 0).any() else 0
        u = np.where(g >= 0, g, G + np.cumsum(g < 0) - 1) + nxt
        nxt = int(u.max()) + 1
        tips += [t['id'] for t in h['tips']]
        unit += u.tolist()
        comp += [rd.name] * len(g)
    log(f'{len(region_dirs)} components, {len(tips)} tips, {nxt} super-units at {UNIT_LEVEL}')
    return tips, np.array(unit), comp


def plan(region_dirs, intake_dir, out, hierarchy_name='hierarchy.json', depth=DEPTH, covers=COVERS, threads=8, log=print, seed=0):
    """Write OUT/{units.json, seqs.fasta, manifest.json} for the given regions. `seed` draws an
    independent plan (other covers, other members) for a stability check."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tips, unit, comp = super_units(region_dirs, hierarchy_name, log)
    U = int(unit.max()) + 1
    members = collections.defaultdict(list)
    for i, u in enumerate(unit.tolist()):
        members[u].append(i)
    corpus = read_fasta(Path(intake_dir) / 'full.derep.fasta')
    size = regions.copies(intake_dir)
    reps = {}
    for u, mm in members.items():
        for i in sorted(mm, key=lambda i: (-size.get(tips[i], 1), tips[i]))[:REPS]:
            reps[f'u{u}|{tips[i]}'] = corpus[tips[i]]
    q, hits = out / 'reps.fasta', out / 'reps_hits.tsv'
    write_fasta(q, reps)
    subprocess.run(['vsearch', '--usearch_global', str(q), '--db', str(q), '--id', str(NBR_ID), '--maxaccepts', str(4 * REPS * SIZE),
                    '--maxrejects', str(4 * REPS * SIZE), '--threads', str(threads), '--userout', f'{hits}.tmp',
                    '--userfields', 'query+target+id', '--quiet'], check=True, stderr=subprocess.DEVNULL)
    sim = collections.defaultdict(dict)
    for line in open(f'{hits}.tmp'):
        a, b, i = line.rstrip('\n').split('\t')
        ua, ub = a.split('|')[0], b.split('|')[0]
        if ua != ub and float(i) > sim[ua].get(ub, 0.0):
            sim[ua][ub] = float(i)
    Path(f'{hits}.tmp').replace(hits)
    ids = [f'u{u}' for u in range(U)]
    nbrs = {u: [v for v, _ in sorted(sim[u].items(), key=lambda x: (-x[1], x[0]))] for u in ids}
    man_u = cover.plan(nbrs, ids, [], covers, depth, SIZE, log=log, seed=seed)
    shards = {}
    for name, s in sorted(man_u['shards'].items()):
        rng = random.Random(f'corpus:{seed}:{name}')
        us = [int(x[1:]) for x in s['members']]
        have = set(s['members'])
        ctx = [int(x[1:]) for x in nbrs.get(s['seed'], []) if x not in have][:CONTEXT]
        shards[name] = {'cover': 0, 'seed': s['seed'], 'units': us,
                        'members': sorted(tips[rng.choice(members[u])] for u in us),
                        'context': sorted(tips[rng.choice(members[u])] for u in ctx)}
    # a cluster of units smaller than a shard gets one shard per cover (each pair COVERS votes):
    # repeat such shards with other random members until each of their units is in COVERS * depth
    target = covers * depth
    held = collections.Counter(u for s in shards.values() for u in s['units'])
    extra = 0
    for name in sorted(shards):
        s = shards[name]
        need = target - min(held[u] for u in s['units'])
        for j in range(max(0, need)):
            rng = random.Random(f'corpus:{seed}:{name}:r{j}')
            ctx_units = [int(x[1:]) for x in nbrs.get(s['seed'], []) if int(x[1:]) not in set(s['units'])][:CONTEXT]
            shards[f'{name}.r{j}'] = {'cover': 0, 'seed': s['seed'], 'units': s['units'], 'repeat_of': name,
                                      'members': sorted(tips[rng.choice(members[u])] for u in s['units']),
                                      'context': sorted(tips[rng.choice(members[u])] for u in ctx_units)}
            extra += 1
        for u in s['units']:
            held[u] += max(0, need)
    if extra:
        log(f'{extra} repeat shards (other random members) for units in clusters smaller than a shard')
    need = sorted({t for s in shards.values() for t in s['members'] + s['context']})
    write_fasta(out / 'seqs.fasta', {t: corpus[t] for t in need}, need)
    json.dump({'unit_level': UNIT_LEVEL, 'tips': tips, 'unit_of_tip': unit.tolist(), 'component_of_tip': comp,
               'components': [Path(r).name for r in region_dirs]}, open(out / 'units.json', 'w'))
    json.dump({'scaffold': [], 'shards': shards, 'corpus': {'units': U, 'covers': covers, 'depth': depth, 'size': SIZE, 'seed': seed,
                                                          'context': CONTEXT, 'neighbour_identity': NBR_ID, 'reps': REPS}},
              open(out / 'manifest.json', 'w'))
    log(f'corpus plan: {U} units, {len(shards)} shards, {len(need)} sequences')
    return len(shards)


def chunk_names(manifest, chunk, of):
    """The shards of chunk `chunk` of `of`: balanced by size squared (L-INS-i cost), largest first."""
    names = sorted(manifest['shards'], key=lambda n: (-len(manifest['shards'][n]['members']) - len(manifest['shards'][n]['context']), n))
    load = [0] * of
    mine = []
    for n in names:
        k = min(range(of), key=lambda i: (load[i], i))
        load[k] += (len(manifest['shards'][n]['members']) + len(manifest['shards'][n]['context'])) ** 2
        if k == chunk:
            mine.append(n)
    return mine


def trees(out, chunk=0, of=1, procs=8, log=print):
    """Build the trees of one chunk of the corpus plan into OUT/trees (cached; twins once)."""
    out = Path(out)
    man = json.load(open(out / 'manifest.json'))
    mine = chunk_names(man, chunk, of)
    sub = dict(man, shards={n: man['shards'][n] for n in mine})
    log(f'chunk {chunk}/{of}: {len(mine)} shards')
    cover.build_trees(sub, read_fasta(out / 'seqs.fasta'), out / 'trees', procs, log=log)


def layer(out, procs=8, log=print):
    """Prepare every planned tree and build the corpus layer (OUT/layer.json)."""
    out = Path(out)
    man = json.load(open(out / 'manifest.json'))
    missing = [n for n in man['shards'] if not (out / 'trees' / f'{n}.nwk.key.json').exists()]
    if missing:
        raise SystemExit(f'{len(missing)} planned trees missing (e.g. {missing[:3]}): build every chunk first')
    shardtrees.export(out, 0, procs, log=log)
    u = json.load(open(out / 'units.json'))
    _, lines, widths = shardtrees.read(out, 0)
    lay = sparse.build_layer(lines, widths, u['tips'], np.array(u['unit_of_tip']), procs=procs, log=log, levels=LEVELS)
    lay.update({'method': METHOD, 'unit_level': UNIT_LEVEL, 'source': f'{len(lines)} corpus shard trees'})
    sparse._write(out / 'layer.json', lay)
    return lay


def release_view(comps, units, layer):
    """What a release needs from the corpus layer. comps: [(tips, {node: [tip ids]}, nodes at the
    component's top (UNIT_LEVEL) level, capped component linkage)] in release order; units, layer:
    units.json and layer.json. Returns a dict:
      levels      the corpus levels
      nodes       corpus groups of 2+ units, one per distinct tip set: tips, kmin, kmax (indices into
                  levels), components, rows {k: evidence}, parent (('C', i) or None)
      continues   {(ci, nd): last corpus level index a component node is still alone at}
      top_parent  {(ci, nd): the group key holding a top-level component node above its last level}
      tip_rows    {tip id: [per corpus level: (group key or None, membership, votes, pull, target, pull votes)]}
                  target: a group key, a tip id (str), or None
      linkage     the corpus dendrogram over all tips ({Z, observed, possible}) and its tip order
    Group keys: (ci, nd) for component nodes, ('C', i) for corpus nodes."""
    levels = [L['level'] for L in layer['levels']]
    t_unit = dict(zip(units['tips'], units['unit_of_tip']))
    U = layer['units']
    unit_tips = collections.defaultdict(list)
    for t, u in t_unit.items():
        unit_tips[u].append(t)
    # each unit: a top-level component node, or a lone tip
    unit_node = {}
    for ci, (tips, tips_of_node, top, _) in enumerate(comps):
        by_set = {frozenset(tips_of_node[nd]): nd for nd in top}
        for t in tips:
            assert t in t_unit, f'tip {t} is in no corpus unit: corpus plan and release differ'
        for u in {t_unit[t] for t in tips}:
            s = frozenset(unit_tips[u])
            assert s <= set(tips), f'unit {u} spans components'
            if s in by_set:
                unit_node[u] = (ci, by_set[s])
            else:
                assert len(s) == 1, f'unit {u} ({len(s)} tips) is no top-level group of its component'
    comp_of_tip = {t: ci for ci, c in enumerate(comps) for t in c[0]}
    nodes, node_of = [], {}
    keys_at = []                                     # per level: per layer group index, its node key
    for k, L in enumerate(layer['levels']):
        keys = []
        for g in L['groups']:
            tt = sorted(t for u in g['units'] for t in unit_tips[u])
            s = frozenset(tt)
            if s not in node_of:
                node_of[s] = len(nodes)
                nodes.append({'tips': tt, 'kmin': k, 'kmax': k, 'components': len({comp_of_tip[t] for t in tt}),
                              'rows': {}, 'parent': None, 'units': g['units']})
            i = node_of[s]
            nodes[i]['kmax'] = k
            nodes[i]['rows'][k] = {x: g[x] for x in ('cohesion', 'votes', 'pull', 'margin', 'held')}
            keys.append(('C', i))
        keys_at.append(keys)
    continues, top_parent = {}, {}
    for u, key in unit_node.items():
        for k, L in enumerate(layer['levels']):
            if L['unit_group'][u] >= 0:
                top_parent[key] = keys_at[k][L['unit_group'][u]]
                break
            continues[key] = k
    for i, nd in enumerate(nodes):
        k = nd['kmax'] + 1
        if k < len(levels):
            L = layer['levels'][k]
            nd['parent'] = keys_at[k][L['unit_group'][nd['units'][0]]]

    def key_of_unit(k, u):
        L = layer['levels'][k]
        if L['unit_group'][u] >= 0:
            return keys_at[k][L['unit_group'][u]]
        return unit_node.get(u)
    tip_rows = {}
    for t, u in t_unit.items():
        rows = []
        for k, L in enumerate(layer['levels']):
            p = L['unit_pull'][u]
            target = None
            if p is not None:
                share, tg, votes, pt = p
                if tg >= 0:
                    target = keys_at[k][tg]
                elif tg == -1:
                    target = unit_node.get(pt) or unit_tips[pt][0]
            rows.append((key_of_unit(k, u), L['unit_membership'][u], L['unit_m_votes'][u],
                         p[0] if p else None, target, p[2] if p else None))
        tip_rows[t] = rows
    return {'levels': levels, 'nodes': nodes, 'continues': continues, 'top_parent': top_parent, 'tip_rows': tip_rows,
            'linkage': corpus_linkage(comps, units, layer)}


def corpus_linkage(comps, units, layer):
    """One dendrogram over every tip: each component's capped linkage (its rows up to UNIT_LEVEL)
    with the corpus layer's over units above. Layer merges below UNIT_LEVEL are drawn just above it.
    Returns ({Z, observed, possible}, tips in index order)."""
    order = [t for c in comps for t in c[0]]
    gi = {t: i for i, t in enumerate(order)}
    n = len(order)
    Z, obs, pos = [], [], []
    root_of_tip = {}
    off = 0
    for tips, _, _, link in comps:
        nc = len(tips)
        new_id = {i: off + i for i in range(nc)}
        parent = {}
        for r, (a, b, ht, sz) in enumerate(link['Z']):
            if ht is None:
                continue                              # above the cap: replaced by the corpus layer
            k = n + len(Z)
            Z.append([new_id[int(a)], new_id[int(b)], ht, sz])
            obs.append(link['observed'][r])
            pos.append(link['possible'][r])
            new_id[nc + r] = k
            parent[new_id[int(a)]] = k
            parent[new_id[int(b)]] = k
        for i in range(nc):
            x = off + i
            while x in parent:
                x = parent[x]
            root_of_tip[tips[i]] = x
        off += nc
    size = collections.Counter()
    hgt = {}
    for t in order:
        size[root_of_tip[t]] += 1
    for k, row in enumerate(Z):
        hgt[n + k] = row[2]
    t_unit = dict(zip(units['tips'], units['unit_of_tip']))
    node = {}
    for t in order:
        u = t_unit[t]
        r = root_of_tip[t]
        assert node.setdefault(u, r) == r, f'unit {u} is not one clade of its component tree'
    floor = UNIT_LEVEL + 1e-6
    U = layer['units']
    for r, (a, b, ht, _) in enumerate(layer['linkage']):
        na, nb = node[int(a)], node[int(b)]
        k = n + len(Z)
        h_ = None if ht is None else max(ht, hgt.get(na, 0.0) or 0.0, hgt.get(nb, 0.0) or 0.0, floor)
        Z.append([na, nb, h_, size[na] + size[nb]])
        obs.append(layer['observed'][r])
        pos.append(layer['possible'][r])
        node[U + r] = k
        size[k] = size[na] + size[nb]
        hgt[k] = h_ if h_ is not None else max(hgt.get(na, 0.0) or 0.0, hgt.get(nb, 0.0) or 0.0)
    return {'Z': Z, 'observed': obs, 'possible': pos}, order
