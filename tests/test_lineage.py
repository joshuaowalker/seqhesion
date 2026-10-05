from seqhesion import lineage


def test_growth_keeps_the_id_and_new_inputs_do_not_count_against_it():
    old = {'A': {'a', 'b', 'c'}}
    new = {'x': {'a', 'b', 'c', 'n1', 'n2', 'n3', 'n4'}}
    inherit, rows = lineage.match(new, old)
    assert inherit == {'x': 'A'}
    assert rows[0]['jaccard'] == 1.0 and rows[0]['event'] == 'continued'


def test_a_split_keeps_the_id_only_for_a_majority_part():
    old = {'A': set('abcdefgh')}
    inherit, _ = lineage.match({'x': set('abcdef'), 'y': set('gh')}, old)
    assert inherit == {'x': 'A'}                       # 6/8 > 1/2
    inherit, rows = lineage.match({'x': set('abcd'), 'y': set('efgh')}, old)
    assert inherit == {}                               # 4/8 is not more than half: both new, A retired
    assert {r['event'] for r in rows} == {'born', 'retired'}


def test_nested_near_duplicates_each_keep_their_own():
    old = {'C': set('abcde'), 'P': set('abcdef')}
    new = {'c': set('abcde'), 'p': set('abcdef')}
    inherit, _ = lineage.match(new, old)
    assert inherit == {'c': 'C', 'p': 'P'}
    inherit, rows = lineage.match({'p': set('abcdef')}, old)   # the child vanished: the parent keeps its id
    assert inherit == {'p': 'P'}
    assert any(r['event'] == 'retired' and r['old'] == 'C' and r['new'] == 'p' for r in rows)


def test_minted_ids_are_reproducible_and_never_reused():
    a = lineage.mint('fp', ['i1', 'i2'], set())
    assert a == lineage.mint('fp', ['i2', 'i1'], set())
    assert len(a) == 9 and a[4] == '-'
    b = lineage.mint('fp', ['i1', 'i2'], {a})
    assert b != a and b.startswith(a.replace('-', '')[:4])


def test_a_parent_adding_only_unshared_inputs_does_not_take_the_childs_id():
    # the coarse parent adds input 'z', which the previous release held in no group: on the shared
    # inputs both new groups equal the old one; the exact one keeps the id
    old = {'OLD': {'a', 'b', 'c'}, 'OTHER': {'q', 'r'}}
    new = {(0, 1): {'a', 'b', 'c'}, (0, 9): {'a', 'b', 'c', 'z'}, (0, 2): {'q', 'r'}}
    inherit, _ = lineage.match(new, old)
    assert inherit[(0, 1)] == 'OLD' and (0, 9) not in inherit
