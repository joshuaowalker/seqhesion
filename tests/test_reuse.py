import json
import os

from seqhesion import cover, tools

SEQS = {'a': 'ACGTACGTAC', 'b': 'ACGTACGTAA', 'c': 'ACGAACGTAC', 'd': 'ACGTTCGTAC'}


def fake_tree(d, name, ids):
    """A complete earlier tree under d: the files of TREE_FILES, keyed by these ids in this order."""
    d.mkdir(parents=True, exist_ok=True)
    for ext in cover.TREE_FILES[:-1]:
        (d / f'{name}{ext}').write_text('(' + ','.join(f'{i}:0.1' for i in ids) + ');' if ext == '.nwk' else 'x')
    json.dump({'key': tools.build_key({i: SEQS[i] for i in ids}, order=list(ids))}, open(d / f'{name}.nwk.key.json', 'w'))


def test_an_earlier_tree_is_reused_by_content_under_the_new_shard_name(tmp_path):
    old = tmp_path / 'old'
    fake_tree(old, 'c0.0007', ['a', 'b', 'c'])            # the same input under another name
    fake_tree(old, 'c0.0000', ['a', 'b', 'd'])            # a different input under the new name
    man = {'scaffold': [], 'shards': {'c0.0000': {'cover': 0, 'seed': 'a', 'members': ['c', 'a', 'b']}}}
    new = tmp_path / 'new'
    assert cover.link_reusable(man, SEQS, new, [old], log=lambda m: None) == 1
    assert os.path.samefile(new / 'c0.0000.nwk', old / 'c0.0007.nwk')
    assert json.load(open(new / 'c0.0000.nwk.key.json'))['key'] == json.load(open(old / 'c0.0007.nwk.key.json'))['key']


def test_nothing_is_linked_without_an_exact_key(tmp_path):
    old = tmp_path / 'old'
    fake_tree(old, 'c0.0000', ['a', 'b'])
    man = {'scaffold': [], 'shards': {'c0.0000': {'cover': 0, 'seed': 'a', 'members': ['a', 'b', 'c']}}}
    assert cover.link_reusable(man, SEQS, tmp_path / 'new', [old], log=lambda m: None) == 0


def test_chunks_partition_the_jobs():
    jobs = [{'region': f'r{k % 3}', 'names': [f'c0.{k:04d}'], 'ids': list('ab' * (k % 5 + 1))} for k in range(17)]
    parts = [cover.chunk_jobs(jobs, c, 4) for c in range(4)]
    names = sorted(j['names'][0] for p in parts for j in p)
    assert names == sorted(j['names'][0] for j in jobs) and all(parts)
