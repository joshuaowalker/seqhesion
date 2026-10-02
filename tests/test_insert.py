from seqhesion import insert

LEVELS = [0.005, 0.01, 0.02]
# two fine groups (a, b) inside one coarse group (C); tip z alone until 0.02
GROUP_AT = {
    0.005: {'a1': 'a', 'a2': 'a', 'b1': 'b', 'b2': 'b', 'z': None},
    0.01: {'a1': 'a', 'a2': 'a', 'b1': 'b', 'b2': 'b', 'z': None},
    0.02: {'a1': 'C', 'a2': 'C', 'b1': 'C', 'b2': 'C', 'z': 'C'},
}


def test_finest_group_within_the_level_then_its_ancestors():
    groups, avg = insert.assign({'a1': 0.003, 'a2': 0.004, 'b1': 0.015, 'b2': 0.016}, GROUP_AT, LEVELS)
    assert groups == ['a', 'a', 'C'] and abs(avg - 0.0035) < 1e-12


def test_too_far_for_fine_groups_joins_only_where_the_level_allows():
    groups, _ = insert.assign({'a1': 0.012, 'a2': 0.013, 'b1': 0.015}, GROUP_AT, LEVELS)
    assert groups == [None, None, 'C']


def test_nothing_within_any_level_is_ungrouped():
    groups, avg = insert.assign({'a1': 0.5, 'b1': 0.6}, GROUP_AT, LEVELS)
    assert groups == [None, None, None] and avg is None


def test_central_shards_prefers_the_seed_then_neighbour_rank():
    man = {'shards': {'s1': {'cover': 0, 'members': ['x', 't']}, 's2': {'cover': 0, 'members': ['t', 'x']},
                      's3': {'cover': 0, 'members': ['y', 'z', 't']}, 's4': {'cover': 1, 'members': ['q', 't']},
                      's5': {'cover': 0, 'members': ['y', 'z']}}}
    assert insert.central_shards(man, 't', per_cover=2) == ['s2', 's1', 's4']
