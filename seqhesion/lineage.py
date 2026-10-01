"""Antenomina: group identities that persist from release to release.

An antenomen ("before the name") is a group of the hierarchy with an identity that outlives a
single release, so that comments and curation can attach to it before (and independently of) any
name. Every distinct group gets one.

Matching a release to the previous one:
  * groups are compared by their INPUTS (mm-to-ref's content hashes), which survive changes to
    our own extraction, and only over inputs present in both releases, so a group that merely
    gains new sequences is not penalised;
  * a new group inherits an old id iff the two are each other's best match (Jaccard) and the
    Jaccard exceeds 1/2 -- most of each is the other;
  * anything else gets a new id, minted once and never reused.
Every comparison is logged (lineage.tsv); what to do with a comment on a group that split or
merged is the curator's decision, not ours.

Ids are 8 Crockford base-32 characters shown as XXXX-XXXX. They look random but are reproducible:
a hash of the release's input fingerprint and the group's inputs, lengthened on a collision with
any id ever minted.
"""
import collections
import hashlib

ALPHABET = '0123456789ABCDEFGHJKMNPQRSTVWXYZ'       # Crockford: no I, L, O, U


def _b32(digest, n):
    x = int.from_bytes(digest, 'big')
    out = []
    for _ in range(n):
        x, r = divmod(x, 32)
        out.append(ALPHABET[r])
    return ''.join(out)


def mint(salt, inputs, taken, n=8):
    """A new id for a group with these inputs, not in `taken`."""
    digest = hashlib.sha256((salt + '|' + ','.join(sorted(inputs))).encode()).digest()
    while True:
        code = _b32(digest, n)
        name = '-'.join(code[i:i + 4] for i in range(0, n, 4))
        if name not in taken:
            return name
        n += 1


def match(new, old):
    """new, old: {key: set of inputs}. Returns (inherit {new key: old key}, rows), where rows hold
    every new group's best old match and every old group's best new match with their overlaps."""
    common = set().union(*new.values()) & set().union(*old.values()) if new and old else set()
    by_input = collections.defaultdict(list)
    for k, s in old.items():
        for i in s & common:
            by_input[i].append(k)
    size_new = {k: len(s & common) for k, s in new.items()}
    size_old = {k: len(s & common) for k, s in old.items()}
    ov = collections.Counter()
    for k, s in new.items():
        for i in s & common:
            for o in by_input[i]:
                ov[(k, o)] += 1
    best_new, best_old = {}, {}
    for (k, o), c in ov.items():
        j = c / (size_new[k] + size_old[o] - c)
        if (j, o) > best_new.get(k, (-1, '')):        # ties: the lexically larger key, deterministically
            best_new[k] = (j, o)
        if (j, k) > best_old.get(o, (-1, '')):
            best_old[o] = (j, k)
    inherit = {k: o for k, (j, o) in best_new.items() if j > 0.5 and best_old[o][1] == k}
    rows = []
    for k, (j, o) in best_new.items():
        c = ov[(k, o)]
        rows.append({'new': k, 'old': o, 'jaccard': j, 'shared': c, 'new_common': size_new[k], 'old_common': size_old[o],
                     'event': 'continued' if inherit.get(k) == o else 'born'})
    for k in new:
        if k not in best_new:
            rows.append({'new': k, 'old': None, 'jaccard': None, 'shared': 0, 'new_common': size_new[k], 'old_common': None,
                         'event': 'born'})
    inherited = set(inherit.values())
    for o in old:
        if o in inherited:
            continue
        if o in best_old:
            j, k = best_old[o]
            rows.append({'new': k, 'old': o, 'jaccard': j, 'shared': ov[(k, o)], 'new_common': size_new[k],
                         'old_common': size_old[o], 'event': 'retired'})
        else:
            rows.append({'new': None, 'old': o, 'jaccard': None, 'shared': 0, 'new_common': None, 'old_common': size_old[o],
                         'event': 'retired'})
    return inherit, rows
