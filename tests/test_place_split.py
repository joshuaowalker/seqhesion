import json

import pytest

from seqhesion import place

LEVELS = [0.01, 0.05]
GROUP_AT = {0.01: {'a': 'G1', 'b': 'G1', 'c': 'G2'}, 0.05: {'a': 'H', 'b': 'H', 'c': 'H'}}


def plan():
    rows = [dict(dict.fromkeys(place.PLACE_COLUMNS), input_id=q) for q in ('q1', 'q2', 'q3')]
    rows[2]['status'] = 'no_match'
    return {'rows': rows, 'its2': {'q1': 'ACGT', 'q2': 'ACGA'},
            'jobs': [['q1', 'C', 's1'], ['q1', 'C', 's2'], ['q2', 'C', 's1']]}


def results():
    return [{'q': 'q1', 'shard': 's1', 'joins': {'a': 0.004, 'b': 0.006, 'c': 0.03}, 'moved': False},
            {'q': 'q1', 'shard': 's2', 'joins': {'a': 0.008, 'c': 0.04}, 'moved': True},
            {'q': 'q2', 'shard': 's1', 'joins': {'a': 0.03, 'b': 0.04, 'c': 0.06}, 'moved': False}]


def test_reduce_pools_medians_and_assigns():
    rows, status = place.reduce(plan(), results(), GROUP_AT, LEVELS, log=lambda m: None)
    by = {r['input_id']: r for r in rows}
    # q1: a = median(0.004, 0.008) = 0.006, b = 0.006 -> G1 average 0.006 <= 0.01
    assert (by['q1']['status'], by['q1']['group_id'], by['q1']['level'], by['q1']['shards'], by['q1']['rearranged']) == \
        ('placed', 'G1', 0.01, 2, 1)
    # q2: G1 average 0.035 > 0.01; at 0.05, H average (0.03 + 0.04 + 0.06) / 3 = 0.0433 <= 0.05
    assert (by['q2']['status'], by['q2']['group_id'], by['q2']['level']) == ('placed', 'H', 0.05)
    assert by['q3']['status'] == 'no_match' and status == {'placed': 2, 'no_match': 1}


def test_reduce_needs_every_result_exactly_once():
    with pytest.raises(SystemExit, match='no result'):
        place.reduce(plan(), results()[:2], GROUP_AT, LEVELS, log=lambda m: None)
    with pytest.raises(SystemExit, match='duplicated'):
        place.reduce(plan(), results() + results()[:1], GROUP_AT, LEVELS, log=lambda m: None)


def test_chunks_partition_the_jobs():
    jobs = [[f'q{k}', 'C', f's{k % 7}'] for k in range(23)]
    parts = [place.chunk_jobs(jobs, c, 5) for c in range(5)]
    assert [j for p in parts for j in p] == jobs and all(parts)


def test_results_must_come_from_one_chunking(tmp_path):
    (tmp_path / 'joins').mkdir()
    for name in ('000-of-002.jsonl', '000-of-003.jsonl'):
        (tmp_path / 'joins' / name).write_text(json.dumps(results()[0]) + '\n')
    with pytest.raises(SystemExit, match='different --of'):
        list(place.read_results(tmp_path))
