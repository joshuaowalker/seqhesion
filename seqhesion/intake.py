"""Intake: raw union → full-ITS sequences, oriented and exactly dereplicated.

  pyitsx delimit   region coordinates, strand, chimera flag (our own step, run on raw input)
  classify         our completeness rule, not pyitsx's full_ITS span (ubertree lab, DESIGN.md §3.1):
                     full   ITS1 >= 50, 5.8S >= 140, ITS2 >= 50, and each flank either found
                            (SSU / LSU anchor) or, if missing, the region beside it >= 150
                            (about the 5th percentile of anchored ITS1/ITS2; Josh 2026-09-22:
                            GenBank submissions are trimmed at the ITS boundaries, so a
                            missing flank there is not truncation). Missing flanks are
                            recorded in classes.tsv `anchors` (SSU+LSU / SSU / LSU / none).
                     its2   LSU anchor, ITS2 >= 100, but ITS1 < 50 or 5.8S < 140
                     other  end-truncated, ITS1-only, undetected
                     chimeric
  extract          full_ITS span of the full class (pyitsx extract orients 5'->3')
  derep            exact dereplication of the extracted sequences; every input id kept

Everything dropped is counted in classes.tsv / summary.json; nothing is corrected.
Usage: python -m seqhesion intake INPUT.fasta OUTDIR
"""
import csv
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict

PYITSX = os.path.join(os.path.dirname(sys.executable), 'pyitsx')
RULE = {'ITS1': 50, '5.8S': 140, 'ITS2': 50}
UNANCHORED_MIN = 150                        # ITS1 (no SSU) or ITS2 (no LSU) must be at least this long
IUPAC = set('ACGTRYSWKMBDHVN')


def span(s):
    if not s or s == '-' or s == 'None':
        return 0
    a, b = s.split('-')
    return int(b) - int(a) + 1


def classify(row):
    if row['chimeric'] == 'True':
        return 'chimeric'
    have = {k: span(row[k]) for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')}
    ok = all(have[k] >= v for k, v in RULE.items())
    ok = ok and (have['SSU'] or have['ITS1'] >= UNANCHORED_MIN) and (have['LSU'] or have['ITS2'] >= UNANCHORED_MIN)
    if ok:
        return 'full'
    if have['LSU'] and have['ITS2'] >= 100:
        return 'its2'
    return 'other'


def anchors(row):
    return '+'.join(k for k in ('SSU', 'LSU') if span(row[k])) or 'none'


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                yield name, ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, ''.join(buf)


def run(union, outdir, cpus=6, log=print):
    os.makedirs(outdir, exist_ok=True)
    readable = os.path.join(outdir, 'readable.fasta')
    unreadable = []
    with open(readable, 'w') as fh:                       # pyitsx rejects anything outside IUPAC
        for name, seq in read_fasta(union):
            if set(seq) - IUPAC:
                unreadable.append(name)
            else:
                fh.write(f'>{name}\n{seq}\n')
    delim = os.path.join(outdir, 'delimit.tsv')
    if not os.path.exists(delim):
        subprocess.run([PYITSX, 'delimit', '-i', readable, '-o', delim + '.tmp', '--cpus', str(cpus)], check=True)
        os.replace(delim + '.tmp', delim)
    rows = list(csv.DictReader(open(delim), delimiter='\t'))
    cls = {r['seq_id']: classify(r) for r in rows}
    counts = Counter(cls.values())
    counts['unreadable'] = len(unreadable)
    counts['full_unanchored'] = sum(1 for r in rows if cls[r['seq_id']] == 'full' and anchors(r) != 'SSU+LSU')
    rows += [{'seq_id': n, 'strand': '', 'chimeric': '', **{k: '' for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')}} for n in unreadable]
    cls.update({n: 'unreadable' for n in unreadable})
    counts['undetected_by_delimit'] = sum(1 for r in rows if not any(r[k] not in ('', '-', 'None') for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')))
    counts['reverse_strand'] = sum(1 for r in rows if r['strand'] == '-')
    with open(os.path.join(outdir, 'classes.tsv'), 'w') as fh:
        fh.write('seq_id\tclass\tstrand\tanchors\tSSU\tITS1\t5.8S\tITS2\tLSU\n')
        for r in rows:
            fh.write('\t'.join([r['seq_id'], cls[r['seq_id']], r['strand'], anchors(r) if cls[r['seq_id']] != 'unreadable' else '']
                               + [str(span(r[k])) for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')]) + '\n')
    log(f'classes: {dict(counts)}')

    full_in = os.path.join(outdir, 'full_input.fasta')
    with open(full_in, 'w') as fh:
        for name, seq in read_fasta(union):
            if cls.get(name) == 'full':
                fh.write(f'>{name}\n{seq}\n')
    ext = os.path.join(outdir, 'full.extract.fasta')
    if not os.path.exists(ext):
        subprocess.run([PYITSX, 'extract', '-i', full_in, '-o', ext + '.tmp', '--cpus', str(cpus), '--region', 'full_ITS'], check=True)
        os.replace(ext + '.tmp', ext)

    groups = defaultdict(list)
    for name, seq in read_fasta(ext):
        groups[seq].append(name.split('|')[0])
    reps = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[1][0]))
    with open(os.path.join(outdir, 'full.derep.fasta'), 'w') as fh, open(os.path.join(outdir, 'derep_members.tsv'), 'w') as mh:
        mh.write('unique\tsize\tmembers\n')
        for k, (seq, names) in enumerate(reps, 1):
            fh.write(f'>u{k:06d}\n{seq}\n')
            mh.write(f'u{k:06d}\t{len(names)}\t{",".join(names)}\n')
    summary = {'input': len(rows), 'classes': dict(counts), 'extracted': sum(len(v) for v in groups.values()),
               'unique': len(groups), 'singletons': sum(len(v) == 1 for v in groups.values()),
               'largest_identical_group': max((len(v) for v in groups.values()), default=0)}
    json.dump(summary, open(os.path.join(outdir, 'summary.json'), 'w'), indent=1)
    log(json.dumps(summary, indent=1))
    return summary


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('union')
    ap.add_argument('outdir')
    ap.add_argument('--cpus', type=int, default=6)
    a = ap.parse_args()
    run(a.union, a.outdir, a.cpus)
