"""An input directory from any FASTA: the form `intake` and `release` read.

  sequences.fasta   one record per distinct sequence, its id the content hash (intake.content_id:
                    min(sha256(s), sha256(revcomp(s))) on the cleaned sequence), sorted by id
  names.tsv         input_id and every header that carried that sequence (display only: nothing
                    in seqhesion reads names)
  manifest.json     schema seqhesion-input/0.1; fingerprint_sha256 = sha256 of sequences.fasta

Exact copies, in either orientation, become one input. Nothing else is changed: sequences are
cleaned (uppercase; whitespace, '-' and '.' removed) and otherwise kept as given.
"""
import datetime
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from .intake import clean, content_id

SCHEMA = 'seqhesion-input/0.1'


def _records(path):
    head, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if head is not None:
                yield head, ''.join(buf)
            head, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if head is not None:
        yield head, ''.join(buf)


def prepare(fasta, out, log=print):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    seqs, names, records, empty = {}, defaultdict(list), 0, 0
    for head, seq in _records(fasta):
        records += 1
        s = clean(seq)
        if not s:
            empty += 1
            continue
        i = content_id(s)
        seqs.setdefault(i, s)
        names[i].append(head)
    with open(out / 'sequences.fasta', 'w') as f:
        for i in sorted(seqs):
            f.write(f'>{i}\n{seqs[i]}\n')
    with open(out / 'names.tsv', 'w') as f:
        f.write('input_id\tname\n')
        for i in sorted(names):
            f.writelines(f'{i}\t{h}\n' for h in names[i])
    counts = {'records': records, 'empty': empty, 'sequences': len(seqs)}
    json.dump({'schema': SCHEMA,
               'created': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
               'fingerprint_sha256': hashlib.sha256((out / 'sequences.fasta').read_bytes()).hexdigest(),
               'input_id': "min(sha256(s), sha256(revcomp(s))) on clean_sequence(s): uppercase, whitespace and "
                           "'-'/'.' stripped; complement ACGTRYSWKMBDHVN->TGCAYRSWMKVHDBN",
               'source': Path(fasta).name, 'counts': counts}, open(out / 'manifest.json', 'w'), indent=1)
    log(f'{records} records -> {len(seqs)} distinct sequences ({empty} empty records skipped)')
    return counts
