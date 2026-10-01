"""Covers that change only locally when sequences are added (experimental, 2026-10-01).

cover.greedy_cover draws seeds from one random stream over the whole region, so adding a few
sequences reshuffles every shard (v0 -> v1: about half of the shards kept their members where 6-11
sequences were added), and with them the trees and ~10% of the group ids. Here:
  * every tip has a fixed pseudo-random priority, sha256(cover:tip); tips are visited in that order
    and a tip still in fewer than `depth` shards seeds one (itself if unused, else its nearest unused
    neighbour whose shard, with the tip put in, is new). Adding a tip changes the neighbourhoods
    near it, and the seeds those touch; the rest of the cover stays as it was;
  * the scaffold is the first n tips in hash order that are not within the neighbour search of
    one already chosen (identity < `far` to every chosen one), so it too changes only if a new
    tip wins a place;
  * a shard is named by its cover and seed, so a shard that did not change keeps its name and its tree.
"""
import hashlib


def priority(cover, tip):
    return hashlib.sha256(f'{cover}:{tip}'.encode()).hexdigest()


def stable_cover(nbrs, all_ids, size, depth, cover, min_size=4):
    """[(seed, members)] in which every tip with at least min_size - 1 neighbours is in >= depth shards
    where its neighbourhood allows."""
    order = sorted(all_ids, key=lambda t: priority(cover, t))
    coverage = dict.fromkeys(all_ids, 0)
    shards, used, seen = [], set(), set()

    def take(seed, members):
        used.add(seed)
        seen.add(frozenset(members))
        shards.append((seed, members))
        for m in members:
            coverage[m] += 1

    for t in order:
        if len(nbrs.get(t, ())) < min_size - 1:
            continue
        tried = 0
        while coverage[t] < depth:
            cands = ([t] if t not in used else []) + [s for s in nbrs[t] if s not in used and len(nbrs.get(s, ())) >= min_size - 1]
            for seed in cands:
                members = [seed] + nbrs[seed][:size - 1]
                if t not in members:
                    members = members[:-1] + [t]
                if frozenset(members) not in seen:
                    take(seed, members)
                    break
            else:
                break                      # no fresh context left: leave it under-covered
            tried += 1
            if tried > depth * 4:
                break
    return shards


def stable_scaffold(all_ids, ident, n=30, far=85.0):
    """n tips in hash order, each below `far` % identity to every one already chosen (pairs absent
    from the neighbour search count as far)."""
    chosen = []
    for t in sorted(all_ids, key=lambda t: priority('scaffold', t)):
        if all(ident.get(t, {}).get(c, 0.0) < far and ident.get(c, {}).get(t, 0.0) < far for c in chosen):
            chosen.append(t)
            if len(chosen) == n:
                break
    return chosen
