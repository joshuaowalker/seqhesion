"""A tree may be reused only if it was made from exactly what is being asked for."""
import pytest

from seqhesion import cache, tools
from seqhesion.fasta import write_fasta

SEQS = {'a': 'ACGTACGTAC', 'b': 'ACGTACGTAA', 'c': 'ACGAACGTAC', 'd': 'ACGTTCGTAC'}


def fake_build(monkeypatch):
    calls = []

    def build(seqs, workdir, name, order=None):
        calls.append(name)
        ids = list(order or seqs)
        write_fasta(workdir / f'{name}.fasta', seqs, ids)
        (workdir / f'{name}.nwk').write_text('(' + ','.join(f'{i}:0.1' for i in ids) + ');')
        return workdir / f'{name}.nwk', {'n': len(ids)}     # same contract as the real build_tree
    monkeypatch.setattr(tools, 'build_tree', build)
    return calls


def test_key_depends_on_content_and_order():
    k = tools.build_key(SEQS)
    assert k == tools.build_key(dict(SEQS))
    assert k != tools.build_key({**SEQS, 'a': 'ACGTACGTAG'})
    assert k != tools.build_key(SEQS, order=['d', 'c', 'b', 'a'])


def test_a_verified_mafft_build_suffix_does_not_enter_the_key(monkeypatch):
    real = cache.tool_version

    def key_with(mafft_version):
        monkeypatch.setattr(cache, 'tool_version', lambda p: mafft_version if p == tools.MAFFT else real(p))
        return tools.build_key(SEQS)
    assert key_with('7.526-opt4') == key_with('7.526')
    assert key_with('7.527') != key_with('7.526')


def test_same_name_different_members_is_refused(tmp_path, monkeypatch):
    calls = fake_build(monkeypatch)
    _, reused = tools.cached_build_tree(SEQS, tmp_path, 'c0.0001')
    assert not reused
    _, reused = tools.cached_build_tree(SEQS, tmp_path, 'c0.0001')
    assert reused and calls == ['c0.0001']
    with pytest.raises(cache.Stale):
        tools.cached_build_tree({**SEQS, 'e': 'ACGTACGTTT'}, tmp_path, 'c0.0001')


def test_unkeyed_tree_is_refused(tmp_path, monkeypatch):
    fake_build(monkeypatch)
    tools.build_tree(SEQS, tmp_path, 'old')
    with pytest.raises(cache.Stale):
        tools.cached_build_tree(SEQS, tmp_path, 'old')


def test_tips_are_checked_even_when_the_key_matches(tmp_path, monkeypatch):
    fake_build(monkeypatch)
    tools.cached_build_tree(SEQS, tmp_path, 't')
    (tmp_path / 't.nwk').write_text('(a:0.1,b:0.1,c:0.1);')
    with pytest.raises(cache.Stale):
        tools.cached_build_tree(SEQS, tmp_path, 't')
