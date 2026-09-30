"""Hierarchy from co-association evidence that covers only some pairs .

A pair of tips is *observed* when at least one shard holds both. In a region of thousands of tips
only a few percent of pairs are, because a shard is a neighbourhood. A pair no shard holds is not
evidence that the two are far apart; it is no evidence at all. So the linkage here averages over
observed pairs only:

  distance(A, B) = mean of d(a, b) over the observed pairs a in A, b in B

and two clusters with no observed pair between them are never merged on evidence. The update is a
weighted mean of the two parents' distances (weights = observed pair counts), so, as in ordinary
average linkage, merge heights never decrease. Each merge records how many of its cross pairs were
observed, so a join made on thin evidence is visible rather than hidden.
"""
import heapq
from multiprocessing import Pool

import numpy as np

UNJOINED = np.inf   # height of the rows that join clusters no observed pair connects


def observed_average_linkage(i, j, d, n):
    """Average linkage over the observed pairs (i[k], j[k]) at distance d[k]; 0 <= i < j < n, each
    pair at most once. Returns (Z, observed, possible):
      Z         scipy-format linkage, n-1 rows; the rows joining parts of a forest (no observed pair
                between them) come last, at height UNJOINED, largest part first
      observed  per row, observed cross pairs behind the merge (0 for UNJOINED rows)
      possible  per row, |A| * |B|
    Ties are broken by the lower slot index, so the result depends only on the input pairs, not
    their order."""
    i, j, d = np.asarray(i, np.int64), np.asarray(j, np.int64), np.asarray(d, np.float64)
    assert np.all(i < j) and (len(i) == 0 or j.max() < n)
    assert len(np.unique(i * n + j)) == len(i), 'each pair at most once'
    adj = [{} for _ in range(n)]
    heap = []
    for a, b, x in zip(i.tolist(), j.tolist(), d.tolist()):
        e = [x, 1]                              # [sum of distances, observed pairs], shared both ways
        adj[a][b] = adj[b][a] = e
        heap.append((x, a, b, 1))
    heapq.heapify(heap)
    sid = list(range(n))                        # slot -> scipy cluster id
    size = [1] * n
    alive = [True] * n
    rows, observed, possible = [], [], []
    while heap:
        x, a, b, c = heapq.heappop(heap)
        if not (alive[a] and alive[b]):
            continue
        e = adj[a].get(b)
        if e is None or e[1] != c:              # stale: the pair has gained evidence since
            continue
        keep, gone = (a, b) if len(adj[a]) >= len(adj[b]) else (b, a)
        rows.append([min(sid[a], sid[b]), max(sid[a], sid[b]), e[0] / e[1], size[a] + size[b]])
        observed.append(e[1])
        possible.append(size[a] * size[b])
        del adj[keep][gone], adj[gone][keep]
        for y, f in adj[gone].items():
            del adj[y][gone]
            g = adj[keep].get(y)
            if g is None:
                adj[keep][y] = adj[y][keep] = f
            else:
                g[0] += f[0]
                g[1] += f[1]
                f = g
            lo, hi = min(keep, y), max(keep, y)
            heapq.heappush(heap, (f[0] / f[1], lo, hi, f[1]))
        adj[gone] = None
        alive[gone] = False
        sid[keep] = n + len(rows) - 1
        size[keep] += size[gone]
    heights = np.array([r[2] for r in rows])
    # a weighted mean of two distances >= the minimum merged; only rounding can dip below
    assert np.all(np.diff(heights) >= -1e-12), 'inversion'
    for k in range(1, len(rows)):
        rows[k][2] = max(rows[k][2], rows[k - 1][2])
    # the forest: parts with no observed pair between them, joined without evidence
    roots = sorted((s for s in range(n) if alive[s]), key=lambda s: (-size[s], sid[s]))
    while len(roots) > 1:
        a, b = roots[0], roots.pop(1)
        rows.append([min(sid[a], sid[b]), max(sid[a], sid[b]), UNJOINED, size[a] + size[b]])
        observed.append(0)
        possible.append(size[a] * size[b])
        sid[a] = n + len(rows) - 1
        size[a] += size[b]
    return np.array(rows, dtype=np.float64).reshape(-1, 4), np.array(observed), np.array(possible)


def forest_roots(Z):
    """Parts of the hierarchy joined by evidence: 1 + the number of UNJOINED rows."""
    return 1 + int(np.sum(Z[:, 2] == UNJOINED))


def votes_at(K, V, level):
    """Per observed pair, from every shard's join level (K = pair keys sorted, V = join levels):
    (unique keys, hits = shards joining the pair at or below `level`, votes = shards holding it)."""
    starts = np.flatnonzero(np.r_[True, K[1:] != K[:-1]])
    hits = np.add.reduceat((V <= level).astype(np.int64), starts)
    votes = np.diff(np.r_[starts, len(K)])
    return K[starts], hits, votes


def pooled(keys, hits, votes, n, group_of):
    """Co-association statistics at one level, POOLED: each shard's verdict on a pair is one vote,
    so a pair seen by twenty shards weighs twenty times one seen by a single shard.

    group_of: per tip, its group (>= 0) or -1 if it is in no group of 2+ (a singleton).
    Returns a dict of arrays:
      membership, m_votes, m_partners   per tip: share of the votes on its pairs with the rest of
                                        its group that join them; votes; partners behind them
      cohesion, c_votes, c_pairs        per group: the same over all its within-group pairs
      pull, p_group, p_tip, p_votes, p_partners
                                        per tip: the outside group (or single ungrouped tip) with
                                        the highest pooled share, ties to more votes then the lower
                                        group; p_group -1 = an ungrouped tip (p_tip), -2 = no
                                        outside evidence at all, -3 = outside evidence, none joining
    Pairs no shard holds cast no votes."""
    G = int(group_of.max()) + 1 if len(group_of) and group_of.max() >= 0 else 0
    i, j = keys // n, keys % n
    gi, gj = group_of[i], group_of[j]
    same = (gi == gj) & (gi >= 0)
    h, v = hits[same], votes[same]
    m_h = np.bincount(i[same], h, n) + np.bincount(j[same], h, n)
    m_v = np.bincount(i[same], v, n) + np.bincount(j[same], v, n)
    m_p = np.bincount(i[same], minlength=n) + np.bincount(j[same], minlength=n)
    c_h = np.bincount(gi[same], h, G)
    c_v = np.bincount(gi[same], v, G)
    c_p = np.bincount(gi[same], minlength=G)
    out = {'membership': np.where(m_v > 0, m_h / np.maximum(m_v, 1), np.nan), 'm_votes': m_v.astype(int),
           'm_partners': m_p, 'cohesion': np.where(c_v > 0, c_h / np.maximum(c_v, 1), np.nan),
           'c_votes': c_v.astype(int), 'c_pairs': c_p}
    # pull: pool each tip's cross-boundary votes by the partner's group (an ungrouped partner is
    # its own target), then keep the best target per tip
    x = ~same
    a, b = np.r_[i[x], j[x]], np.r_[j[x], i[x]]
    hh, vv = np.r_[hits[x], hits[x]], np.r_[votes[x], votes[x]]
    tgt = np.where(group_of[b] >= 0, group_of[b], G + b)
    code = a * (G + n) + tgt
    u, inv = np.unique(code, return_inverse=True)
    H, V, P = np.bincount(inv, hh), np.bincount(inv, vv), np.bincount(inv)
    tip, t = u // (G + n), u % (G + n)
    ratio = H / V
    o = np.lexsort((t, -V, -ratio, tip))
    first = o[np.unique(tip[o], return_index=True)[1]]
    pull = np.zeros(n)
    p_group = np.full(n, -2)
    p_tip = np.full(n, -1)
    p_votes = np.zeros(n, int)
    p_partners = np.zeros(n, int)
    tb, tt = tip[first], t[first]
    pull[tb] = ratio[first]
    p_group[tb] = np.where(tt < G, tt, -1)
    p_tip[tb] = np.where(tt < G, -1, tt - G)
    p_votes[tb] = V[first]
    p_partners[tb] = P[first]
    none = (p_group != -2) & (pull == 0)
    p_group[none] = -3
    # with nothing joining, report all outside evidence, not the arbitrary first target
    if none.any():
        tot_v = np.bincount(a, vv, n).astype(int)
        tot_p = np.bincount(a, minlength=n)
        p_votes[none], p_partners[none], p_tip[none] = tot_v[none], tot_p[none], -1
    out.update(pull=pull, p_group=p_group, p_tip=p_tip, p_votes=p_votes, p_partners=p_partners)
    return out


VERDICTS = ('clade', 'unresolved', 'conflict', 'whole')


def join_matrix(tree, index, join='diameter', floor=0.0):
    """(member tip indices, matrix of per-pair join levels) for one shard tree. 'diameter': the
    diameter of the pair's MRCA clade; 'compatible': at a polytomy, the diameter of just the two
    children holding them; 'identical': that only when it is <= floor (the data's resolution), else
    the diameter (see joins below)."""
    height, diam, below = {}, {}, {}
    leaves = [l for l in tree.iter_leaves()]
    pos = {l: k for k, l in enumerate(leaves)}
    J = np.zeros((len(leaves), len(leaves)))
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            height[node], diam[node], below[node] = 0.0, 0.0, [pos[node]]
            continue
        reach = sorted((height[c] + c.dist for c in node.children), reverse=True)
        height[node] = reach[0]
        diam[node] = max(max(diam[c] for c in node.children), reach[0] + (reach[1] if len(reach) > 1 else 0.0))
        kids = node.children
        for a in range(len(kids)):
            for b in range(a + 1, len(kids)):
                ca, cb = kids[a], kids[b]
                v = diam[node] if join == 'diameter' else max(diam[ca], diam[cb], height[ca] + ca.dist + height[cb] + cb.dist)
                if join == 'identical' and v > floor:
                    v = diam[node]
                J[np.ix_(below[ca], below[cb])] = v
                J[np.ix_(below[cb], below[ca])] = v
        below[node] = [k for c in kids for k in below.pop(c)]
    return [index[l.name] for l in leaves], J, height, diam


def shard_evidence(tree, index, width, groups_of_tip, group_mask, group_level, levels, join='diameter'):
    """What one shard tree says about each group it holds >= 2 tips of.

    tree           an ete3 tree of the shard's members (rooted, scaffold pruned, as exported)
    groups_of_tip  per tip index, the groups (ids) containing it
    group_mask     per group, its tips as a bitmask
    group_level    per group, the index into `levels` at which to count joins (its display level)
    join           the per-pair join level the hierarchy was built with ('diameter', 'compatible')
    Returns [(group, held, outsiders, verdict, stem, hits, votes)]:
      verdict   'clade' the held tips form a clade of this tree; 'conflict' a clade of the tree
                takes some of them together with an outsider; 'unresolved' neither (a polytomy
                holds them with outsiders); 'whole' the shard holds no outsider, so it cannot tell
      stem      for 'clade': the branch above it in expected changes (length x alignment width)
      hits, votes   pairs of held tips whose join level in this tree is <= the group's level, of
                all held pairs -- this shard's contribution to cohesion"""
    tip_ids, J, _, _ = join_matrix(tree, index, join, 1.0 / width)
    where = {t: k for k, t in enumerate(tip_ids)}
    mask, clades = {}, {}
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            mask[node] = 1 << index[node.name]
            continue
        m = 0
        for c in node.children:
            m |= mask[c]
        mask[node] = m
        if not node.is_root():
            clades[m] = node.dist
    members = mask[tree]
    count = {}
    for t in tip_ids:
        for g in groups_of_tip[t]:
            count[g] = count.get(g, 0) + 1
    n_members = len(tip_ids)
    rows = []
    for g, held in count.items():
        if held < 2:
            continue
        gm = group_mask[g]
        h = gm & members
        out = n_members - held
        stem = None
        if out == 0:
            verdict = 'whole'
        elif h in clades:
            verdict, stem = 'clade', clades[h] * width
        elif any((c & h) and (c & h) != h and (c & ~gm & members) for c in clades):
            verdict = 'conflict'
        else:
            verdict = 'unresolved'
        pos = [where[t] for t in tip_ids if gm >> t & 1]
        sub = J[np.ix_(pos, pos)][np.triu_indices(len(pos), 1)]
        hits = int(np.sum(sub <= levels[group_level[g]]))
        rows.append((g, held, out, verdict, stem, hits, held * (held - 1) // 2))
    return rows


def separation(keys, dist, n, tips, group_of):
    """Spread and nearest relatives of one group, from consensus pairwise distances.

    keys, dist  observed pairs (i*n+j, i<j) and their distance (median over shards)
    tips        the group's tips
    group_of    per tip, its group at the level used to name relatives (-1 = ungrouped)
    Returns {'spread', 'spread_median', 'within_pairs',
             'nearest': (distance, tip, group) -- the nearest outside GROUP of 2+ tips, via its
                        closest tip; the "closest sibling", which may sit anywhere in the hierarchy,
             'second':  the same for the next group, or None,
             'stray':   (distance, tip) the nearest ungrouped outside tip, or None,
             'strays_closer': ungrouped tips nearer than the nearest group}.
    Only tips outside the group are relatives, so its own subgroups never are. Ungrouped tips are
    reported apart: in real regions they are often sequences near-identical to members that
    collapsed polytomies left unjoined (39-54% of groups had one as their nearest relative).
    Pairs no shard holds are not distances and are ignored."""
    tips = np.asarray(tips)
    inside = np.zeros(n, bool)
    inside[tips] = True
    i, j = keys // n, keys % n
    ii, jj = inside[i], inside[j]
    w = dist[ii & jj]
    x = ii ^ jj
    other = np.where(ii[x], j[x], i[x])
    d = dist[x]
    res = {'spread': float(w.max()) if len(w) else None, 'spread_median': float(np.median(w)) if len(w) else None,
           'within_pairs': int(len(w)), 'nearest': None, 'second': None, 'stray': None, 'strays_closer': 0}
    if not len(d):
        return res
    g = group_of[other]
    o = np.lexsort((other, d))
    grouped = o[g[o] >= 0]
    if len(grouped):
        first = grouped[np.unique(g[grouped], return_index=True)[1]]
        first = first[np.lexsort((other[first], d[first]))][:2]
        rel = [(float(d[k]), int(other[k]), int(g[k])) for k in first]
        res['nearest'] = rel[0]
        res['second'] = rel[1] if len(rel) > 1 else None
    loose = o[g[o] < 0]
    if len(loose):
        res['stray'] = (float(d[loose[0]]), int(other[loose[0]]))
        limit = res['nearest'][0] if res['nearest'] else np.inf
        res['strays_closer'] = int(len(np.unique(other[loose][d[loose] < limit])))
    return res


# ---------------------------------------------------------------------------------------------
# Per-shard join levels and their consensus over a cover.

_G = {}


# 'identical' is the adopted join (2026-09-30); 'diameter' and 'compatible' are kept for comparison,
# 'patristic' feeds the distance annotations.
JOIN_VALUES = ('diameter', 'compatible', 'identical', 'height', 'patristic')


def _init(index, value='diameter'):
    _G.update(index=index, value=value)


def _joins(job):
    line, width = job
    _G['floor'] = 1.0 / width if width else 0.0
    return joins(line)


def joins(line):
    """(keys i*n+j with i<j, join level) for every pair of member tips in one exported tree.
    'diameter': the diameter of their most recent common ancestor's clade. 'height': twice the
    MRCA's height (its distance down to its farthest tip), so the level is on about the diameter
    scale. A parent's height exceeds its child's by at least the child's stem, while its diameter
    need not exceed the child's at all (a clade wider than its stem): height keeps the stem.
    'compatible': the tightest group holding both that the tree does not contradict: at a
    resolved node the same as 'diameter'; at a polytomy, the diameter of just the two children
    holding them, since any union of a polytomy's children is compatible with it (a polytomy is
    missing information, not conflict -- Page 2002). With 'diameter' every pair in a polytomy
    joins at the polytomy's full width, which split near-identical tips between groups
    (ubertree lab, experiments/coassoc/RESULTS.md, 2026-09-30).
    'identical': like 'compatible', but only for pairs whose span through the polytomy is below
    the data's resolution (one expected change = 1 / alignment width): indistinguishable tips join
    at their span, every other pair in a polytomy at its full width. 'compatible' groups the
    distinguishable lineages inside a polytomy by distance, which the simulation showed is mostly
    noise at that scale (lab simulation: wrong groups 4.4% -> 8.6% at 0.03).
    'patristic': not a join level but the pair's own tip-to-tip distance in this tree (for the
    nearest-relative and spread annotations; median over shards like the others)."""
    from ete3 import Tree
    t = Tree(line, format=5)
    index, n = _G['index'], len(_G['index'])
    height, diam, leaves, ldepth = {}, {}, {}, {}
    depth = {}
    for node in t.traverse('preorder'):
        depth[node] = depth[node.up] + node.dist if node.up is not None else 0.0
    ks, vs = [], []
    value = _G['value']
    for node in t.traverse('postorder'):
        if node.is_leaf():
            height[node], diam[node] = 0.0, 0.0
            leaves[node] = np.array([index[node.name]], dtype=np.int64)
            ldepth[node] = np.array([depth[node]])
            continue
        reach = sorted((height[c] + c.dist for c in node.children), reverse=True)
        height[node] = reach[0]
        diam[node] = max(max(diam[c] for c in node.children), reach[0] + (reach[1] if len(reach) > 1 else 0.0))
        kids = [leaves[c] for c in node.children]
        kd = [ldepth[c] for c in node.children]
        for a in range(len(kids)):
            for b in range(a + 1, len(kids)):
                x, y = np.meshgrid(kids[a], kids[b], indexing='ij')
                lo, hi = np.minimum(x, y).ravel(), np.maximum(x, y).ravel()
                ks.append(lo * n + hi)
                if value == 'patristic':
                    vs.append((kd[a][:, None] + kd[b][None, :] - 2 * depth[node]).ravel().astype(np.float32))
                elif value in ('compatible', 'identical'):
                    ca, cb = node.children[a], node.children[b]
                    v = max(diam[ca], diam[cb], height[ca] + ca.dist + height[cb] + cb.dist)
                    if value == 'identical' and v > _G['floor']:
                        v = diam[node]
                    vs.append(np.full(lo.size, v, dtype=np.float32))
                else:
                    vs.append(np.full(lo.size, diam[node] if value == 'diameter' else 2 * height[node], dtype=np.float32))
        leaves[node] = np.concatenate(kids)
        ldepth[node] = np.concatenate(kd)
        for c in node.children:
            del leaves[c], ldepth[c]
    return np.concatenate(ks), np.concatenate(vs)


def consensus_levels(lines, widths, index, procs, value='identical'):
    """(keys, median join level, shards holding the pair, (K, V) every shard's value) over one
    cover's prepared shard trees (`shardtrees.read`). `widths` (trimmed alignment columns) set
    the resolution floor for the 'identical' join."""
    assert value in JOIN_VALUES
    K, V = [], []
    with Pool(procs, initializer=_init, initargs=(index, value)) as pool:
        for k, v in pool.imap_unordered(_joins, list(zip(lines, widths)), chunksize=16):
            K.append(k)
            V.append(v)
    K, V = np.concatenate(K), np.concatenate(V)
    order = np.lexsort((V, K))
    K, V = K[order], V[order]
    starts = np.flatnonzero(np.r_[True, K[1:] != K[:-1]])
    counts = np.diff(np.r_[starts, len(K)])
    med = (V[starts + (counts - 1) // 2] + V[starts + counts // 2]) / 2
    return K[starts], med, counts, (K, V)


def hierarchy(keys, level, n):
    """Average-linkage hierarchy on the consensus join levels, over observed pairs only ->
    (Z, observed, possible)."""
    return observed_average_linkage(keys // n, keys % n, level, n)
