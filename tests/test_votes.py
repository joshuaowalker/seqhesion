import numpy as np

from seqhesion import corpus, export
from seqhesion.coassoc import observed_average_linkage


def test_merges_carry_the_votes_on_their_observed_cross_pairs():
    # 0-1 close on 5 votes; 2 joins them through two pairs on 1 and 2 votes; 3 is seen by nothing
    Z, obs, pos, votes = observed_average_linkage([0, 0, 1], [1, 2, 2], [0.1, 0.2, 0.3], 4, w=[5, 1, 2])
    assert Z[0, :2].tolist() == [0, 1] and votes[0] == 5
    assert obs[1] == 2 and pos[1] == 2 and votes[1] == 3
    assert votes[2] == 0                                      # the forest join: no evidence


def test_without_weights_the_result_is_unchanged():
    a = observed_average_linkage([0, 0, 1], [1, 2, 2], [0.1, 0.2, 0.3], 3)
    b = observed_average_linkage([0, 0, 1], [1, 2, 2], [0.1, 0.2, 0.3], 3, w=[1, 1, 1])
    assert len(a) == 3 and np.array_equal(a[0], b[0]) and b[3].tolist() == b[1].tolist()


def test_a_group_reports_its_join_and_its_weakest_merge():
    # tips a b c d: a+b on 9 votes, +c on 1 vote, +d on 20 votes (layer rows); a fine row has no votes
    Z = [[0, 1, 0.1, 2], [2, 4, 0.2, 3], [3, 5, 0.4, 4]]
    rows = [{'group_id': 'G1', 'kind': 'level', 'parent_id': 'G2'}, {'group_id': 'G2', 'kind': 'level', 'parent_id': 'G3'},
            {'group_id': 'G3', 'kind': 'level', 'parent_id': None}]
    tips_of = {'G1': ['a', 'b'], 'G2': ['a', 'b', 'c'], 'G3': ['a', 'b', 'c', 'd']}
    export.node_columns(rows, tips_of, [({'Z': Z, 'votes': [None, 1, 20]}, ['a', 'b', 'c', 'd'])])
    by = {r['group_id']: r for r in rows}
    assert (by['G1']['join_votes'], by['G1']['weakest_votes']) == (None, None)       # fine merge only
    assert (by['G2']['join_votes'], by['G2']['weakest_votes']) == (1, 1)
    assert (by['G3']['join_votes'], by['G3']['weakest_votes']) == (20, 1)            # a thin merge inside


def test_a_verification_round_is_chunked_on_its_own():
    man = {'shards': {'c0.0001': {'members': list('abc'), 'context': []},
                      'v1.0000': {'members': list('ab'), 'context': [], 'verify_round': 1},
                      'v1.0001': {'members': list('abcd'), 'context': [], 'verify_round': 1},
                      'v2.0000': {'members': list('ab'), 'context': [], 'verify_round': 2}}}
    got = corpus.chunk_names(man, 0, 2, round_=1) + corpus.chunk_names(man, 1, 2, round_=1)
    assert sorted(got) == ['v1.0000', 'v1.0001']
    assert len(corpus.chunk_names(man, 0, 1)) == 4
