from seqhesion import intake


def row(**kw):
    r = {'chimeric': 'False', 'SSU': '1-44', 'ITS1': '45-294', '5.8S': '295-452', 'ITS2': '453-650', 'LSU': '651-688'}
    r.update(kw)
    return r


def test_completeness_rule():
    assert intake.classify(row()) == 'full'
    assert intake.classify(row(chimeric='True')) == 'chimeric'
    assert intake.classify(row(ITS1='45-55', **{'5.8S': '56-100'})) == 'its2'      # gITS7-style: ITS2 present, ITS1/5.8S stubs
    assert intake.classify(row(SSU='')) == 'full'                                    # flankless but ITS1 250: complete
    assert intake.classify(row(SSU='', ITS1='45-150')) == 'its2'                     # no SSU and ITS1 only 106: truncated
    assert intake.classify(row(LSU='-')) == 'full'                                   # no LSU, ITS2 198
    assert intake.classify(row(LSU='-', ITS2='453-580')) == 'other'                  # no LSU, ITS2 128
    assert intake.classify(row(SSU='', LSU='')) == 'full' and intake.anchors(row(SSU='', LSU='')) == 'none'
    assert intake.anchors(row()) == 'SSU+LSU' and intake.anchors(row(LSU='')) == 'SSU'


def test_content_id_is_orientation_free_and_iupac_aware():
    s = 'ACGTRYSWKMBDHVN'
    rc = s.translate(intake.COMP)[::-1]
    assert rc == 'NBDHVKMWSRYACGT'
    assert intake.content_id(s) == intake.content_id(rc)
    assert intake.content_id(intake.clean(' acg-t.\nRY ')) == intake.content_id('ACGTRY')


def test_input_ids_must_be_their_content_hash(tmp_path):
    import pytest
    good = intake.content_id('ACGTACGGT')
    p = tmp_path / 'in.fasta'
    p.write_text(f'>{good} label="anything at all"\nacgt-acggt\n')
    assert intake.checked_input(p) == [(good, 'ACGTACGGT')]
    p.write_text(f'>{good}\nACGTACGGA\n')
    with pytest.raises(intake.BadInput):
        intake.checked_input(p)
    p.write_text(f'>{good}\nACGTACGGT\n>{good}\nACGTACGGT\n')
    with pytest.raises(intake.BadInput):
        intake.checked_input(p)


def test_tip_id_is_the_extracted_content():
    assert intake.tip_id('acgt') == intake.tip_id('ACGT') != intake.tip_id('ACGA')
    assert len(intake.tip_id('ACGT')) == intake.TIP_HEX
