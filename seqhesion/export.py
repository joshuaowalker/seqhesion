"""A release: the hierarchy of one or more components as TSV tables plus a JSON manifest.

Schema seqhesion-release/0 (agreed with mm-to-ref, 2026-09-30). No labels anywhere.

  manifest.json    schema, seqhesion commit, the input (directory + fingerprint), method, levels,
                   components built, counts, and the columns of every table
  inputs.tsv       input_id, tip_id, status, component
                   status: built (its tip is in a component this release built) | not_built |
                   dropped:<intake class> (its2, other, chimeric, unreadable, overlong, not_extracted)
  tips.fasta       tip_id and its extracted, oriented full-ITS sequence, for every tip in a built component
  tips.tsv         tip_id, component, n_inputs, seen (held by some shard), part (forest part, 0 = largest)
  groups.tsv       one row per distinct group (the same tips over a run of levels), measured at its
                   coarsest level: group_id, component, n_tips, level_min, level_max, parent_id,
                   cohesion, votes, replication, votes_replicate, held, stem, stem_shards,
                   spread, spread_median, identity_min, identity_median, nearest_id, nearest_distance,
                   second_id, second_distance, gap, clade, unresolved, conflict, whole
  group_levels.tsv group_id, level, cohesion, votes, replication, votes_replicate, held
  group_members.tsv group_id, tip_id
  lineage.tsv      event (continued / born / retired), group_id, previous_id, jaccard, shared,
                   new_common, old_common: every group's best match in the previous release and
                   every previous group's best match here, over the inputs both releases hold
  placements.tsv   every ITS2-only input (intake class its2), placed through the shard trees
                   (seqhesion.place, seqhesion.insert): input_id, status (placed / ungrouped /
                   not_built / spans_components / no_match / no_its2), component, group_id and level
                   (the finest group it joins; its ancestors follow from groups.tsv), join (its
                   average join level to that group, tree distance), shards (shard trees it was
                   placed in), rearranged (of those, how many re-inferences moved its surroundings: a
                   strain signal), candidate and candidate_identity (the tip and vsearch ITS2
                   identity used only to choose where to look). Placements never shape the hierarchy.
  minted.tsv       every group id ever minted, with the release that minted it (never reused)
  membership.tsv   tip_id, level, group_id (empty: in no group), membership, votes, partners,
                   pull, pull_target (a group_id, or a tip_id for an ungrouped tip; empty when nothing
                   outside joins it at this level), pull_votes; pull empty: no outside evidence at all

Tree data (schema 0.1, for drawing and on-demand distances; see the release manifest's notes):
  dendrogram/<component>.nwk          the average-linkage dendrogram the levels are cut from (primary
                                      cover): one Newick line per forest part, internal nodes m<row>,
                                      branch lengths = height differences (ultrametric)
  dendrogram/<component>.merges.tsv   row, height, size, observed, possible (cross pairs behind it)
  shards/<component>/cover<c>.nwk.gz  the prepared shard trees as the hierarchy reads them, one per line
  shards/<component>/cover<c>.shards.tsv  shard, members held, branches collapsed, trimmed width

group_id is the group's antenomen: an identity carried from release to release (seqhesion.lineage:
a group inherits the id of its mutual best match in the previous release when they share more
than half their inputs; otherwise a new id is minted). component is the tip_id of the component's first tip (most inputs,
then tip_id). Levels, spreads and distances are tree distances (substitutions per site, medians
over shard trees); identity is vsearch global pairwise identity (--iddef 2), a different scale.
"""
import collections
import csv
import datetime
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from . import lineage, place
from .fasta import read_fasta, write_fasta

SCHEMA = 'seqhesion-release/0.1'
COLUMNS = {
    'inputs.tsv': ['input_id', 'tip_id', 'status', 'component'],
    'tips.tsv': ['tip_id', 'component', 'n_inputs', 'seen', 'part'],
    'groups.tsv': ['group_id', 'component', 'n_tips', 'level_min', 'level_max', 'parent_id', 'cohesion', 'votes',
                   'replication', 'votes_replicate', 'held', 'stem', 'stem_shards', 'spread', 'spread_median',
                   'identity_min', 'identity_median', 'nearest_id', 'nearest_distance', 'second_id', 'second_distance',
                   'gap', 'clade', 'unresolved', 'conflict', 'whole'],
    'group_levels.tsv': ['group_id', 'level', 'cohesion', 'votes', 'replication', 'votes_replicate', 'held'],
    'group_members.tsv': ['group_id', 'tip_id'],
    'lineage.tsv': ['event', 'group_id', 'previous_id', 'jaccard', 'shared', 'new_common', 'old_common'],
    'minted.tsv': ['group_id', 'release'],
    'placements.tsv': place.PLACE_COLUMNS,
    'membership.tsv': ['tip_id', 'level', 'group_id', 'membership', 'votes', 'partners', 'pull', 'pull_target', 'pull_votes'],
}


def _fmt(x):
    if x is None:
        return ''
    if isinstance(x, float):
        return f'{x:.6g}'
    return str(x)


class Table:
    def __init__(self, path, name):
        self.f = open(path / name, 'w')
        self.cols = COLUMNS[name]
        self.f.write('\t'.join(self.cols) + '\n')

    def row(self, **kw):
        assert set(kw) == set(self.cols), set(kw) ^ set(self.cols)
        self.f.write('\t'.join(_fmt(kw[c]) for c in self.cols) + '\n')

    def close(self):
        self.f.close()


def pairwise_identity(seqs, threads=8):
    """{(a, b): identity in [0, 1]} for every pair of `seqs` ({id: sequence}), vsearch global."""
    if len(seqs) < 2:
        return {}
    with tempfile.TemporaryDirectory() as d:
        write_fasta(f'{d}/in.fasta', seqs)
        subprocess.run(['vsearch', '--allpairs_global', f'{d}/in.fasta', '--acceptall', '--iddef', '2', '--threads', str(threads),
                        '--userout', f'{d}/out.tsv', '--userfields', 'query+target+id', '--quiet'],
                       stderr=subprocess.DEVNULL, check=True)
        out = {}
        for line in open(f'{d}/out.tsv'):
            q, t, i = line.rstrip('\n').split('\t')
            out[(q, t) if q < t else (t, q)] = float(i) / 100
    return out


def identities(h, tips, seqs, threads=8, log=print):
    """Per node: (min, median) pairwise identity over its tips. Pairs are aligned within each
    group of the coarsest level only; every node lies inside one of those."""
    top = h['levels'][-1]
    node_of_top = {}
    res = {}
    for g in top['groups']:
        ids = [tips[t] for t in g['tips']]
        ident = pairwise_identity({i: seqs[i] for i in ids}, threads)
        node_of_top[g['node']] = (ids, ident)
    # every node: find its top-level group through any of its tips
    top_of_tip = {t: g['node'] for g in top['groups'] for t in g['tips']}
    tips_of_node = collections.defaultdict(list)
    for lv in h['levels']:
        for g in lv['groups']:
            tips_of_node[g['node']] = g['tips']
    for node, tt in tips_of_node.items():
        ids, ident = node_of_top[top_of_tip[tt[0]]]
        names = sorted(tips[t] for t in tt)
        v = [ident[(a, b)] for i, a in enumerate(names) for b in names[i + 1:]]
        res[node] = (min(v), float(np.median(v))) if v else (None, None)
    log(f'identity: {len(res)} groups, {sum(len(x[1]) for x in node_of_top.values()):,} pairs aligned')
    return res


def seqhesion_commit():
    try:
        root = Path(__file__).resolve().parents[1]
        rev = subprocess.run(['git', '-C', str(root), 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', str(root), 'status', '--porcelain', '--', 'seqhesion'], capture_output=True, text=True).stdout.strip()
        return rev + ('-dirty' if dirty else '')
    except (OSError, subprocess.CalledProcessError):
        return 'unknown'


def read_tsv(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


def previous_groups(prev):
    """{group_id: set of input_ids} and the ids ever minted, from a previous release directory."""
    prev = Path(prev)
    inputs_of = collections.defaultdict(set)
    for r in read_tsv(prev / 'inputs.tsv'):
        if r['tip_id']:
            inputs_of[r['tip_id']].add(r['input_id'])
    groups = collections.defaultdict(set)
    for r in read_tsv(prev / 'group_members.tsv'):
        groups[r['group_id']] |= inputs_of[r['tip_id']]
    minted = [(r['group_id'], r['release']) for r in read_tsv(prev / 'minted.tsv')]
    return dict(groups), minted


def dendrogram(linkage, tips):
    """(Newick lines, one per forest part; merge rows) from hierarchy.build's linkage."""
    Z, obs, pos = linkage['Z'], linkage['observed'], linkage['possible']
    n = len(tips)
    height = {i: 0.0 for i in range(n)}
    children = {}
    for k, (a, b, h, size) in enumerate(Z):
        children[n + k] = (a, b)
        height[n + k] = h
    roots = []
    for k, (a, b, h, size) in enumerate(Z):        # forest parts: rows at height None join them
        if h is None:
            for c in (a, b):
                if height[c] is not None:
                    roots.append(c)
    if not roots:
        roots = [n + len(Z) - 1] if Z else list(range(n))
    seen_roots = []
    for r in roots:                                  # a part may itself be a None-joined row's child only once
        if r not in seen_roots:
            seen_roots.append(r)
    lines = []
    for r in sorted(seen_roots, key=lambda r: -(Z[r - n][3] if r >= n else 1)):
        out, stack = [], [(r, 'open')]
        while stack:                                 # iterative: the tree can be deep
            node, state = stack.pop()
            parent_h = None
            if state == 'open':
                if node < n:
                    out.append(('leaf', node))
                else:
                    a, b = children[node]
                    stack.append((node, 'close'))
                    stack.append((b, 'open'))
                    stack.append((node, 'comma'))
                    stack.append((a, 'open'))
                    out.append(('(', node))
            elif state == 'comma':
                out.append((',', node))
            else:
                out.append((')', node))
        parent_of = {}
        for node, (a, b) in children.items():
            parent_of[a] = node
            parent_of[b] = node

        def blen(node):
            p = parent_of.get(node)
            if p is None or height[p] is None:
                return 0.0
            return height[p] - height[node]
        txt = []
        for kind, node in out:
            if kind == 'leaf':
                txt.append(f'{tips[node]}:{blen(node):.6g}')
            elif kind == '(':
                txt.append('(')
            elif kind == ',':
                txt.append(',')
            else:
                lab = f'm{node - n}'
                txt.append(f'){lab}' + (f':{blen(node):.6g}' if node != r else ''))
        lines.append(''.join(txt) + ';')
    rows = [{'row': k, 'height': h, 'size': size, 'observed': obs[k], 'possible': pos[k]} for k, (a, b, h, size) in enumerate(Z)]
    return lines, rows


def write_tree_data(out, comp, region_dir, h, tips, covers='covers'):
    """dendrogram/ and shards/ for one component (schema 0.1)."""
    import gzip
    import shutil
    from . import shardtrees
    dd = out / 'dendrogram'
    dd.mkdir(exist_ok=True)
    if 'linkage' in h:
        lines, rows = dendrogram(h['linkage'], tips)
        (dd / f'{comp}.nwk').write_text('\n'.join(lines) + '\n')
        with open(dd / f'{comp}.merges.tsv', 'w') as f:
            f.write('row\theight\tsize\tobserved\tpossible\n')
            for r in rows:
                f.write(f"{r['row']}\t{_fmt(r['height'])}\t{r['size']}\t{r['observed']}\t{r['possible']}\n")
    sd = out / 'shards' / comp
    sd.mkdir(parents=True, exist_ok=True)
    cd = Path(region_dir) / covers
    for c in (0, 1):
        nwk, idx = shardtrees.paths(cd, c)
        with open(nwk, 'rb') as src, gzip.open(sd / f'cover{c}.nwk.gz', 'wb') as dst:
            shutil.copyfileobj(src, dst)
        names, _, widths = shardtrees.read(cd, c)
        held = [l.rstrip('\n').split('\t') for l in open(idx)]
        with open(sd / f'cover{c}.shards.tsv', 'w') as f:
            f.write('shard\tmembers_held\tcollapsed\twidth\n')
            for (name, m, col), w in zip(held, widths):
                f.write(f'{name}\t{m}\t{col}\t{w}\n')


def next_name(releases_dir, kind, now=None):
    """The next release name: YYYYMMDD.NN plus 'f' (full build) or 'i' (increment), UTC date, NN the
    release's number that day (Josh, 2026-10-02; provenance lives in the manifest, not the name)."""
    assert kind in ('f', 'i')
    now = now or datetime.datetime.now(datetime.timezone.utc)
    day = now.strftime('%Y%m%d')
    taken = [p.name for p in Path(releases_dir).glob(f'{day}.*')] if Path(releases_dir).exists() else []
    n = 1 + max((int(t[9:11]) for t in taken if len(t) >= 11 and t[9:11].isdigit()), default=0)
    return f'{day}.{n:02d}{kind}'


def write_release(out, input_dir, intake_dir, components, all_components, threads=8, log=print, settings=None, previous=None,
                  place_its2=True, covers='covers'):
    """components: [(region_dir, hierarchy dict from hierarchy.build or its JSON)] to publish.
    all_components: every component of the intake (lists of tip ids), for inputs.tsv.
    previous: the previous release directory, whose group ids are carried forward."""
    out = Path(out)
    release = out.name
    input_dir, intake_dir = Path(input_dir), Path(intake_dir)
    inp = json.load(open(input_dir / 'manifest.json'))
    classes = read_tsv(intake_dir / 'classes.tsv')
    inputs_of_tip = collections.defaultdict(set)
    for r in classes:
        if r['tip']:
            inputs_of_tip[r['tip']].add(r['seq_id'])
    derep = read_fasta(intake_dir / 'full.derep.fasta')
    n_inputs = {t: len(v) for t, v in inputs_of_tip.items()}
    comp_of = {}
    for ids in all_components:
        anchor = min(ids, key=lambda t: (-n_inputs[t], t))
        for t in ids:
            comp_of[t] = anchor

    # pass 1: every group of every component, by its tips
    built_comps = []
    for region_dir, h in components:
        tips = sorted(read_fasta(Path(region_dir) / 'comp.fasta'))
        assert [t['id'] for t in h['tips']] == tips, f'{region_dir}: hierarchy tips differ from the region'
        comp = comp_of[tips[0]]
        assert all(comp_of[t] == comp for t in tips), f'{region_dir} is not one component of this intake'
        tips_of_node = {}
        for lv in h['levels']:
            for g in lv['groups']:
                tips_of_node[g['node']] = [tips[t] for t in g['tips']]
        built_comps.append((region_dir, h, tips, comp, tips_of_node))

    # antenomina: carry ids forward from the previous release, mint the rest
    new = {(ci, nd): set().union(*(inputs_of_tip[t] for t in tt)) for ci, (_, _, _, _, tn) in enumerate(built_comps) for nd, tt in tn.items()}
    old, minted = previous_groups(previous) if previous else ({}, [])
    inherit, rows = lineage.match(new, old)
    taken = {g for g, _ in minted}
    gid = {}
    # deterministic minting order: groups form a hierarchy (any two are nested or disjoint), so
    # (smallest input, size) tells them apart
    for key in sorted(new, key=lambda k: (min(new[k]), len(new[k]))):
        if key in inherit:
            gid[key] = inherit[key]
        else:
            gid[key] = lineage.mint(inp.get('fingerprint_sha256', ''), new[key], taken)
            taken.add(gid[key])
            minted.append((gid[key], release))
    events = collections.Counter(r['event'] for r in rows)
    log(f'antenomina: {len(inherit)} carried forward, {len(new) - len(inherit)} new; events {dict(events)}')

    out.mkdir(parents=True, exist_ok=False)            # a release is written once
    tables = {name: Table(out, name) for name in COLUMNS if name not in ('inputs.tsv', 'lineage.tsv', 'minted.tsv', 'placements.tsv')}
    built_group_at = collections.defaultdict(dict)
    tip_seqs, built = {}, set()
    meta = {'levels': None, 'components': []}
    for ci, (region_dir, h, tips, comp, tips_of_node) in enumerate(built_comps):
        assert meta['levels'] in (None, [lv['level'] for lv in h['levels']])
        meta['levels'] = [lv['level'] for lv in h['levels']]
        built.update(tips)
        tip_seqs.update({t: derep[t] for t in tips})
        node_at = [{g['id']: g['node'] for g in lv['groups']} for lv in h['levels']]
        G = lambda nd: gid[(ci, nd)]  # noqa: E731
        ident = identities(h, tips, derep, threads, log)
        for t in h['tips']:
            tables['tips.tsv'].row(tip_id=t['id'], component=comp, n_inputs=n_inputs[t['id']], seen=int(t['seen']), part=t['part'])
        for nd, tt in tips_of_node.items():
            for t in tt:
                tables['group_members.tsv'].row(group_id=G(nd), tip_id=t)
        for nd_i, nd in enumerate(h['nodes']):
            k = nd['kmax']
            lv = h['levels'][k]
            g = next(x for x in lv['groups'] if x['node'] == nd_i)
            parent = node_at[k + 1][g['parent']] if g['parent'] is not None else None
            near, sec, v = nd['nearest'], nd['second'], nd['verdicts']
            tables['groups.tsv'].row(
                group_id=G(nd_i), component=comp, n_tips=len(tips_of_node[nd_i]),
                level_min=h['levels'][nd['kmin']]['level'], level_max=lv['level'],
                parent_id=G(parent) if parent is not None else None,
                cohesion=g['cohesion'], votes=g['votes'], replication=g['replication'], votes_replicate=g['votes1'],
                held=g['held'], stem=nd['stem'], stem_shards=nd['stem_n'], spread=nd['spread'], spread_median=nd['spread_median'],
                identity_min=ident[nd_i][0], identity_median=ident[nd_i][1],
                nearest_id=G(node_at[k][near[2]]) if near else None, nearest_distance=near[0] if near else None,
                second_id=G(node_at[k][sec[2]]) if sec else None, second_distance=sec[0] if sec else None,
                gap=nd['gap'], clade=v[0], unresolved=v[1], conflict=v[2], whole=v[3])
        for k, lv in enumerate(h['levels']):
            for i, t in enumerate(tips):
                gi = lv['group_of'][i]
                built_group_at[lv['level']][t] = G(node_at[k][gi]) if gi >= 0 else None
            for g in lv['groups']:
                tables['group_levels.tsv'].row(group_id=G(g['node']), level=lv['level'], cohesion=g['cohesion'], votes=g['votes'],
                                               replication=g['replication'], votes_replicate=g['votes1'], held=g['held'])
            for i, t in enumerate(tips):
                gi = lv['group_of'][i]
                p = lv['pull'][i]
                target = None
                if p is not None:
                    if p[1] >= 0:
                        target = G(node_at[k][p[1]])
                    elif p[1] == -1:
                        target = tips[p[4]]
                tables['membership.tsv'].row(
                    tip_id=t, level=lv['level'], group_id=G(node_at[k][gi]) if gi >= 0 else None,
                    membership=lv['membership'][i], votes=lv['m_votes'][i], partners=lv['m_partners'][i],
                    pull=p[0] if p is not None else None, pull_target=target,
                    pull_votes=p[2] if p is not None else None)
        write_tree_data(out, comp, region_dir, h, tips, covers)
        meta['components'].append({'component': comp, 'tips': len(tips), 'groups': len(h['nodes']), 'region_dir': str(region_dir),
                                   'forest': h['forest'], 'provenance': h['provenance']})
        meta['method'], meta['join'] = h['method'], h['join']
        log(f'component {comp}: {len(tips)} tips, {len(h["nodes"])} groups')
    for t in tables.values():
        t.close()
    write_fasta(out / 'tips.fasta', tip_seqs, sorted(tip_seqs))
    lin = Table(out, 'lineage.tsv')
    for r in sorted(rows, key=lambda r: (r['event'], str(r['new']), str(r['old']))):
        lin.row(event=r['event'], group_id=gid[r['new']] if r['new'] is not None else None, previous_id=r['old'],
                jaccard=r['jaccard'], shared=r['shared'], new_common=r['new_common'], old_common=r['old_common'])
    lin.close()
    mt = Table(out, 'minted.tsv')
    for g, rel in minted:
        mt.row(group_id=g, release=rel)
    mt.close()

    placed = None
    if place_its2:
        covers_of = {c['component']: Path(c['region_dir']) / covers for c in meta['components']}
        prow, placed = place.place_inputs(intake_dir, classes, built_group_at, meta['levels'], comp_of, covers_of, threads, log)
        pt = Table(out, 'placements.tsv')
        for r in prow:
            pt.row(**r)
        pt.close()

    status = collections.Counter()
    inputs = Table(out, 'inputs.tsv')
    for r in classes:
        t = r['tip']
        s = ('built' if t in built else 'not_built') if r['class'] == 'full' else f'dropped:{r["class"]}'
        status[s] += 1
        inputs.row(input_id=r['seq_id'], tip_id=t or None, status=s, component=comp_of.get(t) if t else None)
    inputs.close()

    manifest = {
        'schema': SCHEMA,
        'release': release,
        'kind': 'full',                         # every release so far is a full build; increments come next
        'base_full_build': release,
        'increments_since_full': 0,
        'previous_release': str(previous) if previous else None,
        'created': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'seqhesion_commit': seqhesion_commit(),
        'input': {'dir': str(input_dir), 'schema': inp.get('schema'), 'fingerprint_sha256': inp.get('fingerprint_sha256')},
        'method': meta.get('method'), 'join': meta.get('join'), 'levels': meta['levels'], 'settings': settings or {},
        'units': {'level, spread, distance': 'tree distance, substitutions per site (median over shard trees)',
                  'stem': 'expected changes (branch length x trimmed alignment columns), median over shards',
                  'identity': 'vsearch --allpairs_global --iddef 2, fraction'},
        'tree_data': {'dendrogram': 'average linkage over observed pairs of the per-pair median join level (identical join; '
                                    'clade-diameter scale, not branch length); draw ultrametric; structure below ~0.005 is weak',
                      'shards': 'prepared shard trees (rooted on the scaffold, scaffold pruned, branches under half an expected '
                                'change collapsed); cover 0 is the primary; patristic = FastTree GTR+gamma path length on '
                                'trimmed columns; one tree is one sample, the hierarchy uses medians over shards',
                      'stability': 'drawings change with every full rebuild (a new cover); antenomina carry across'},
        'group_id': 'antenomen: carried from the previous release on a mutual best match sharing > 1/2 of the inputs both '
                    'releases hold, else newly minted (see lineage.tsv, minted.tsv)',
        'placements': placed,
        'lineage': {'carried': len(inherit), 'new': len(new) - len(inherit), 'events': dict(events)},
        'components': meta['components'],
        'counts': {'inputs': sum(status.values()), 'status': dict(status), 'tips': len(tip_seqs),
                   'groups': sum(c['groups'] for c in meta['components'])},
        'tables': COLUMNS,
    }
    json.dump(manifest, open(out / 'manifest.json', 'w'), indent=1)
    log(f'wrote release {out}: {dict(status)}')
    return manifest
