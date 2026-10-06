import json

from seqhesion import inputs, intake


def test_prepare_dereplicates_both_orientations(tmp_path):
    fa = tmp_path / 'in.fasta'
    fa.write_text('>a first\nACGT-TTGA\n>b\nacgtttga\n>c\nTCAAACGT\n>d\nGGGCCC\n>e\n\n')
    counts = inputs.prepare(fa, tmp_path / 'out', log=lambda m: None)
    assert counts == {'records': 5, 'empty': 1, 'sequences': 2}
    ids = [l[1:].strip() for l in open(tmp_path / 'out' / 'sequences.fasta') if l.startswith('>')]
    assert ids == sorted(ids) and intake.content_id('ACGTTTGA') in ids
    names = (tmp_path / 'out' / 'names.tsv').read_text().splitlines()
    assert sum(1 for l in names if l.startswith(intake.content_id('ACGTTTGA'))) == 3     # a, b and c (its revcomp)
    assert [r[0] for r in intake.checked_input(tmp_path / 'out' / 'sequences.fasta')] == ids
    assert json.load(open(tmp_path / 'out' / 'manifest.json'))['schema'] == 'seqhesion-input/0.1'
