"""A release: the hierarchy of one or more components as TSV tables plus a JSON manifest.

Schema seqhesion-release/0 (agreed with mm-to-ref, 2026-09-30); 0.2 (2026-10-04): the hierarchy
pools both covers, and replication / votes_replicate are replaced by pull_max / margin; 0.3 (same
day): pull_max (a maximum over members, which saturated in large groups) is replaced by pull, the
members' outside pulls pooled over votes, and margin = cohesion - pull; 0.4 (2026-10-05): with
release --layer, coarse levels 0.125 ... 0.5 from the coarse layer (seqhesion.sparse) follow the
fine ones in every table. No labels anywhere.

Coarse levels (> 0.1): the units are the fine groups at 0.1 (and lone tips). A coarse group of
2+ units is a new group: cohesion / votes / pull / margin / held count UNIT-PAIR votes from the
sparse shards (held = share of its unit pairs observed); stem, spread, identity, nearest and the
verdict counts are empty (no tree-level evidence for it). A fine group still alone at a coarse
level is the same group continuing: its groups.tsv row keeps its fine-level evidence, level_max
reaches its coarsest level, and its group_levels.tsv rows there carry only pull (no unit pairs
inside it). membership.tsv at a coarse level gives each tip its unit's values (partners empty).
The dendrogram is the fine linkage inside each unit and the layer's above them, so every group at
every level is a clade of it; a layer merge below 0.1 is drawn just above it.

0.5 (2026-10-05, agreed with mm-to-ref): with release --corpus, CORPUS levels (0.4, 0.5, 0.75;
manifest corpus_levels) relate components (seqhesion.corpus): their units are each component's
groups at 0.3 (its coarse layer stops there) and lone tips. A corpus group of 2+ units is new; its
`component` is empty when it spans several (n_components, a new groups.tsv column, says how many),
its evidence counts unit-pair votes as for coarse groups. A component's top group still alone at
a corpus level continues (same id, level_max beyond 0.3), and parent_id runs on through the corpus
levels. dendrogram/corpus.nwk (+ .merges.tsv) is one tree of every built tip: the component trees
up to 0.3 and the corpus layer above, every group at every level a clade of it; dendrogram/
<component>.nwk stops at 0.3 (a forest of its 0.3 groups).

0.6 (2026-10-06, agreed with mm-to-ref, adopted by Josh): STEM NODES. A linkage node at height <=
0.1 that is a group at no level (it lives between two levels) but whose linkage stem is >= 0.001
(seqhesion.hierarchy.stem_nodes) is a group too: kind 'node', level_min / level_max empty, cohesion
/ votes / pull / margin / held measured at measured_at (a 0.0005 grid point inside its life), every
other evidence column empty; no group_levels.tsv or membership.tsv rows. Every groups.tsv row gains
kind (level | node), height (its dendrogram node's merge height), height_top (its parent's; empty at
a root), linkage_stem (height_top - height), measured_at (nodes only) and dendrogram_node (its label,
m<row>, in dendrogram/corpus.nwk when the release has corpus levels, else in dendrogram/<component>.nwk).
parent_id is the nearest enclosing group OR node. Stem-node ids carry across releases by the same
lineage rule as every group.

  manifest.json    schema, seqhesion commit, the input (directory + fingerprint), method, levels,
                   components built, counts, and the columns of every table
  inputs.tsv       input_id, tip_id, status, component
                   status: built (its tip is in a component this release built) | not_built |
                   dropped:<intake class> (its2, other, chimeric, unreadable, overlong, not_extracted)
  tips.fasta       tip_id and its extracted, oriented full-ITS sequence, for every tip in a built component
  tips.tsv         tip_id, component, n_inputs, seen (held by some shard), part (forest part, 0 = largest)
  groups.tsv       one row per distinct group (the same tips over a run of levels), measured at its
                   coarsest level: group_id, component, n_tips, level_min, level_max, parent_id,
                   cohesion, votes, pull, margin, held, stem, stem_shards,
                   spread, spread_median, identity_min, identity_median, nearest_id, nearest_distance,
                   second_id, second_distance, gap, clade, unresolved, conflict, whole
  group_levels.tsv group_id, level, cohesion, votes, pull, margin, held
                   cohesion: share of the shard votes on the group's pairs that join them at or
                   below the level; pull: over every member, the votes joining it to its strongest
                   outside target over the votes on those pairs (pooled, like cohesion); margin =
                   cohesion - pull, the confidence measure (seqhesion.hierarchy)
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
  dendrogram/<component>.nwk          the average-linkage dendrogram the levels are cut from: one
                                      Newick line per forest part, internal nodes m<row>,
                                      branch lengths = height differences (ultrametric)
  dendrogram/<component>.merges.tsv   row, height, size, observed, possible (cross pairs behind it)
  shards/<component>/cover<c>.nwk.gz  the prepared shard trees as the hierarchy reads them, one per line
  shards/<component>/cover<c>.shards.tsv  shard, members held, branches collapsed, trimmed width,
                                      twin_of (cover 1: the cover-0 shard with the same members, so the
                                      same tree, counted once by the hierarchy; else empty)

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

SCHEMA = 'seqhesion-release/0.6'
COLUMNS = {
    'inputs.tsv': ['input_id', 'tip_id', 'status', 'component'],
    'tips.tsv': ['tip_id', 'component', 'n_inputs', 'seen', 'part'],
    'groups.tsv': ['group_id', 'component', 'n_components', 'n_tips', 'level_min', 'level_max', 'parent_id', 'cohesion', 'votes',
                   'pull', 'margin', 'held', 'stem', 'stem_shards', 'spread', 'spread_median',
                   'identity_min', 'identity_median', 'nearest_id', 'nearest_distance', 'second_id', 'second_distance',
                   'gap', 'clade', 'unresolved', 'conflict', 'whole',
                   'kind', 'height', 'height_top', 'linkage_stem', 'measured_at', 'dendrogram_node'],
    'group_levels.tsv': ['group_id', 'level', 'cohesion', 'votes', 'pull', 'margin', 'held'],
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


NODE_COLUMNS = ('kind', 'height', 'height_top', 'linkage_stem', 'measured_at', 'dendrogram_node')


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


class GroupRows:
    """groups.tsv, held back until every group is known: kind, heights, dendrogram node and parent_id
    (nearest enclosing group or node) come from the dendrogram (`node_columns`)."""
    def __init__(self, path):
        self.path, self.rows, self.tips = path, [], {}

    def row(self, tips, kind='level', measured_at=None, **kw):
        assert set(kw) == set(COLUMNS['groups.tsv']) - set(NODE_COLUMNS), set(kw) ^ (set(COLUMNS['groups.tsv']) - set(NODE_COLUMNS))
        self.rows.append(dict(kw, kind=kind, measured_at=measured_at))
        self.tips[kw['group_id']] = tips

    def close(self):
        pass

    def write(self, trees):
        node_columns(self.rows, self.tips, trees)
        t = Table(self.path, 'groups.tsv')
        for r in self.rows:
            t.row(**r)
        t.close()


def node_columns(rows, tips_of, trees):
    """Fill each group row's dendrogram columns and parent_id. trees: [(linkage dict, tip order)], each
    tip in exactly one; every group must be a node (clade) of its tree. A level group's parent_id was
    its group at the next coarser level; where a stem node now lies between, it becomes that node."""
    where, info = {}, []
    for ti, (link, order) in enumerate(trees):
        n = len(order)
        Z = link['Z']
        parent = [None] * (n + len(Z))
        size = [1] * n + [int(r[3]) for r in Z]
        height = [0.0] * n + [r[2] for r in Z]
        for r, row in enumerate(Z):
            parent[int(row[0])] = parent[int(row[1])] = n + r
        info.append((n, parent, size, height))
        for i, t in enumerate(order):
            where[t] = (ti, i)
    item = {}
    for r in rows:
        tips = tips_of[r['group_id']]
        ti, x = where[tips[0]]
        n, parent, size, height = info[ti]
        while size[x] < len(tips):
            x = parent[x]
        assert size[x] == len(tips) and x >= n, (r['group_id'], 'is not a node of the dendrogram')
        item[(ti, x)] = r['group_id']
        r['_node'] = (ti, x)
    kinds = {r['group_id']: r['kind'] for r in rows}
    for r in rows:
        ti, x = r['_node']
        n, parent, size, height = info[ti]
        p = parent[x]
        r['dendrogram_node'] = f'm{x - n}'
        r['height'] = height[x]
        r['height_top'] = height[p] if p is not None else None
        r['linkage_stem'] = (height[p] - height[x]) if p is not None and height[p] is not None and height[x] is not None else None
        while p is not None and (ti, p) not in item:
            p = parent[p]
        new = item[(ti, p)] if p is not None else None
        if r['kind'] == 'level' and new != r['parent_id']:
            assert kinds.get(new) == 'node', (r['group_id'], 'parent', r['parent_id'], 'became', new)
        r['parent_id'] = new
        del r['_node']


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
    group of the coarsest FINE level only (a coarse level above it would mean all pairs of a whole
    component); every fine node lies inside one of those, and a node born in the coarse layer
    (sparse.add_layer) gets (None, None)."""
    top = [lv for lv in h['levels'] if not lv.get('layer')][-1]
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
        if lv.get('layer'):
            continue
        for g in lv['groups']:
            tips_of_node[g['node']] = g['tips']
    res.update({k: (None, None) for k, nd in enumerate(h['nodes']) if nd.get('layer')})
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


def write_corpus(out, cx, tables, gid, built_comps, built_group_at, log=print):
    """The corpus levels' rows (groups, members, group levels, membership) and dendrogram/corpus.*."""
    comp_of_tip = {t: comp for _, _, tips, comp, _ in built_comps for t in tips}
    for i, nd in enumerate(cx['nodes']):
        key = ('C', i)
        r = nd['rows'][nd['kmax']]
        comps = {comp_of_tip[t] for t in nd['tips']}
        tables['groups.tsv'].row(
            tips=nd['tips'],
            group_id=gid[key], component=next(iter(comps)) if len(comps) == 1 else None, n_components=len(comps),
            n_tips=len(nd['tips']), level_min=cx['levels'][nd['kmin']], level_max=cx['levels'][nd['kmax']],
            parent_id=gid[nd['parent']] if nd['parent'] else None,
            cohesion=r['cohesion'], votes=r['votes'], pull=r['pull'], margin=r['margin'], held=r['held'],
            stem=None, stem_shards=None, spread=None, spread_median=None, identity_min=None, identity_median=None,
            nearest_id=None, nearest_distance=None, second_id=None, second_distance=None, gap=None,
            clade=None, unresolved=None, conflict=None, whole=None)
        for t in nd['tips']:
            tables['group_members.tsv'].row(group_id=gid[key], tip_id=t)
        for k, rr in sorted(nd['rows'].items()):
            tables['group_levels.tsv'].row(group_id=gid[key], level=cx['levels'][k], **rr)
    for t, rows in cx['tip_rows'].items():
        for k, (key, m, v, pull, target, pv) in enumerate(rows):
            L = cx['levels'][k]
            built_group_at[L][t] = gid[key] if key else None
            tables['membership.tsv'].row(
                tip_id=t, level=L, group_id=gid[key] if key else None, membership=m, votes=v, partners=None, pull=pull,
                pull_target=(gid[target] if isinstance(target, tuple) else target), pull_votes=pv)
    (link, order) = cx['linkage']
    lines, rows = dendrogram(link, order)
    dd = out / 'dendrogram'
    dd.mkdir(exist_ok=True)
    (dd / 'corpus.nwk').write_text('\n'.join(lines) + '\n')
    with open(dd / 'corpus.merges.tsv', 'w') as f:
        f.write('row\theight\tsize\tobserved\tpossible\n')
        for r in rows:
            f.write(f"{r['row']}\t{_fmt(r['height'])}\t{r['size']}\t{r['observed']}\t{r['possible']}\n")
    log(f'corpus: {len(cx["nodes"])} groups, dendrogram/corpus.nwk {len(lines)} part(s) over {len(order)} tips')


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
    from .cover import twins
    twin_of = twins(json.load(open(cd / 'manifest.json')))
    for c in (0, 1):
        nwk, idx = shardtrees.paths(cd, c)
        with open(nwk, 'rb') as src, gzip.open(sd / f'cover{c}.nwk.gz', 'wb') as dst:
            shutil.copyfileobj(src, dst)
        names, _, widths = shardtrees.read(cd, c)
        held = [l.rstrip('\n').split('\t') for l in open(idx)]
        with open(sd / f'cover{c}.shards.tsv', 'w') as f:
            f.write('shard\tmembers_held\tcollapsed\twidth\ttwin_of\n')
            for (name, m, col), w in zip(held, widths):
                f.write(f'{name}\t{m}\t{col}\t{w}\t{twin_of.get(name, "")}\n')


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
                  place_its2=True, covers='covers', corpus=None):
    """components: [(region_dir, hierarchy dict from hierarchy.build or its JSON)] to publish.
    all_components: every component of the intake (lists of tip ids), for inputs.tsv.
    previous: the previous release directory, whose group ids are carried forward.
    corpus: (units.json, layer.json) of the corpus layer (seqhesion.corpus), whose levels follow
    the components' own (their layers capped at its unit level)."""
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

    cx = None
    if corpus:
        from . import corpus as corpus_mod
        cx = corpus_mod.release_view([(tips, tn, {g['node'] for g in h['levels'][-1]['groups']}, h['linkage'])
                                      for _, h, tips, _, tn in built_comps], *corpus)
        log(f'corpus levels {cx["levels"]}: {len(cx["nodes"])} corpus groups')

    # antenomina: carry ids forward from the previous release, mint the rest
    new = {(ci, nd): set().union(*(inputs_of_tip[t] for t in tt)) for ci, (_, _, _, _, tn) in enumerate(built_comps) for nd, tt in tn.items()}
    if cx:
        new.update({('C', i): set().union(*(inputs_of_tip[t] for t in nd['tips'])) for i, nd in enumerate(cx['nodes'])})
    # stem nodes (schema 0.6): groups that live between two levels
    stems = {(ci, 'S', j): [tips[i] for i in x['tips']]
             for ci, (_, h, tips, _, _) in enumerate(built_comps) for j, x in enumerate(h.get('stem_nodes', []))}
    new.update({key: set().union(*(inputs_of_tip[t] for t in tt)) for key, tt in stems.items()})
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
    tables = {name: Table(out, name) for name in COLUMNS
              if name not in ('inputs.tsv', 'lineage.tsv', 'minted.tsv', 'placements.tsv', 'groups.tsv')}
    tables['groups.tsv'] = GroupRows(out)
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
            # measured at its coarsest fine level (kmax); a fine group still alone in the coarse
            # layer continues up to kmax_top, where its parent is
            k = nd['kmax']
            kt = nd.get('kmax_top', k)
            lv = h['levels'][k]
            g = next(x for x in lv['groups'] if x['node'] == nd_i)
            gt = next(x for x in h['levels'][kt]['groups'] if x['node'] == nd_i)
            parent = (ci, node_at[kt + 1][gt['parent']]) if gt['parent'] is not None else None
            level_max = h['levels'][kt]['level']
            if cx and gt['parent'] is None:              # top of its component: the corpus levels go on
                if (ci, nd_i) in cx['continues']:
                    level_max = cx['levels'][cx['continues'][(ci, nd_i)]]
                parent = cx['top_parent'].get((ci, nd_i))
            near, sec, v = nd['nearest'], nd['second'], nd['verdicts'] or [None] * 4
            tables['groups.tsv'].row(
                tips=tips_of_node[nd_i],
                group_id=G(nd_i), component=comp, n_components=1, n_tips=len(tips_of_node[nd_i]),
                level_min=h['levels'][nd['kmin']]['level'], level_max=level_max,
                parent_id=gid[parent] if parent is not None else None,
                cohesion=g['cohesion'], votes=g['votes'], pull=g['pull'], margin=g['margin'],
                held=g['held'], stem=nd['stem'], stem_shards=nd['stem_n'], spread=nd['spread'], spread_median=nd['spread_median'],
                identity_min=ident[nd_i][0], identity_median=ident[nd_i][1],
                nearest_id=G(node_at[k][near[2]]) if near else None, nearest_distance=near[0] if near else None,
                second_id=G(node_at[k][sec[2]]) if sec else None, second_distance=sec[0] if sec else None,
                gap=nd['gap'], clade=v[0], unresolved=v[1], conflict=v[2], whole=v[3])
        for j, x in enumerate(h.get('stem_nodes', [])):
            tt = stems[(ci, 'S', j)]
            tables['groups.tsv'].row(
                tips=tt, kind='node', measured_at=x['measured_at'],
                group_id=gid[(ci, 'S', j)], component=comp, n_components=1, n_tips=len(tt), level_min=None, level_max=None,
                parent_id=None, cohesion=x['cohesion'], votes=x['votes'], pull=x['pull'], margin=x['margin'], held=x['held'],
                stem=None, stem_shards=None, spread=None, spread_median=None, identity_min=None, identity_median=None,
                nearest_id=None, nearest_distance=None, second_id=None, second_distance=None, gap=None,
                clade=None, unresolved=None, conflict=None, whole=None)
            for t in tt:
                tables['group_members.tsv'].row(group_id=gid[(ci, 'S', j)], tip_id=t)
        for k, lv in enumerate(h['levels']):
            for i, t in enumerate(tips):
                gi = lv['group_of'][i]
                built_group_at[lv['level']][t] = G(node_at[k][gi]) if gi >= 0 else None
            for g in lv['groups']:
                tables['group_levels.tsv'].row(group_id=G(g['node']), level=lv['level'], cohesion=g['cohesion'], votes=g['votes'],
                                               pull=g['pull'], margin=g['margin'], held=g['held'])
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
        if cx:
            for key, kc in cx['continues'].items():
                if key[0] != ci:
                    continue
                t0 = tips_of_node[key[1]][0]
                for k in range(kc + 1):
                    tables['group_levels.tsv'].row(group_id=gid[key], level=cx['levels'][k], cohesion=None, votes=0,
                                                   pull=cx['tip_rows'][t0][k][3], margin=None, held=None)
        meta['components'].append({'component': comp, 'tips': len(tips), 'groups': len(h['nodes']),
                                   'stem_nodes': len(h.get('stem_nodes', [])), 'region_dir': str(region_dir),
                                   'forest': h['forest'], 'provenance': h['provenance']})
        meta['method'], meta['join'] = h['method'], h['join']
        if h.get('layer'):
            meta['layer'] = h['layer']['method']
        log(f'component {comp}: {len(tips)} tips, {len(h["nodes"])} groups, {len(h.get("stem_nodes", []))} stem nodes')
    if cx:
        write_corpus(out, cx, tables, gid, built_comps, built_group_at, log)
        meta['levels'] = meta['levels'] + cx['levels']
        meta['corpus_levels'] = cx['levels']
    for t in tables.values():
        t.close()
    tables['groups.tsv'].write([cx['linkage']] if cx else [(h['linkage'], tips) for _, h, tips, _, _ in built_comps])
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
    else:
        Table(out, 'placements.tsv').close()      # the table is always there; empty when placement was skipped

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
                                'change collapsed); the hierarchy pools both covers, a cover-1 shard with twin_of set counted once; patristic = FastTree GTR+gamma path length on '
                                'trimmed columns; one tree is one sample, the hierarchy uses medians over shards',
                      'stability': 'drawings change with every full rebuild (a new cover); antenomina carry across'},
        'group_id': 'antenomen: carried from the previous release on a mutual best match sharing > 1/2 of the inputs both '
                    'releases hold, else newly minted (see lineage.tsv, minted.tsv)',
        'placements': placed if place_its2 else 'skipped (ITS2-only inputs stay dropped:its2)',
        'coarse_layer': meta.get('layer'),
        'corpus_levels': meta.get('corpus_levels'),
        'corpus_layer': corpus[1]['method'] if corpus else None,
        'lineage': {'carried': len(inherit), 'new': len(new) - len(inherit), 'events': dict(events)},
        'components': meta['components'],
        'counts': {'inputs': sum(status.values()), 'status': dict(status), 'tips': len(tip_seqs),
                   'groups': len(tables['groups.tsv'].rows),
                   'stem_nodes': sum(1 for r in tables['groups.tsv'].rows if r['kind'] == 'node')},
        'stem_nodes': ('groups of kind node: linkage nodes at height <= 0.1 that are a group at no level, with linkage '
                       'stem (parent height - own height) >= 0.001; evidence at measured_at, a 0.0005 grid point inside '
                       'their life; dendrogram_node labels are in dendrogram/corpus.nwk when the release has corpus levels, '
                       'else in dendrogram/<component>.nwk'),
        'tables': COLUMNS,
    }
    json.dump(manifest, open(out / 'manifest.json', 'w'), indent=1)
    log(f'wrote release {out}: {dict(status)}')
    return manifest
