"""Intake: input sequences → full-ITS tips, oriented and exactly dereplicated.

Input: a FASTA whose ids are content hashes (the contract with mm-to-ref, 2026-09-30):
  input_id = min(sha256(s), sha256(revcomp(s))), lowercase hex, where s is the sequence
  uppercased with whitespace and gap characters ('-', '.') removed, and revcomp is
  IUPAC-aware. Every id is recomputed here and a mismatch or duplicate refuses the whole
  input. Anything after the id on a header line is ignored (labels are never read).
Tips: tip_id = the first 16 hex digits of sha256 of the extracted, oriented full-ITS
  sequence, so a tip's id is its content (collisions are checked, never assumed away).

  pyitsx delimit   region coordinates, strand, chimera flag (our own step, run on raw input)
  classify         our completeness rule, not pyitsx's full_ITS span (set in the lab prototype):
                     full   ITS1 >= 50, 5.8S >= 140, ITS2 >= 50, and each flank either found
                            (SSU / LSU anchor) or, if missing, the region beside it >= 150
                            (about the 5th percentile of anchored ITS1/ITS2; Josh 2026-09-22:
                            GenBank submissions are trimmed at the ITS boundaries, so a
                            missing flank there is not truncation). Missing flanks are
                            recorded in classes.tsv `anchors` (SSU+LSU / SSU / LSU / none).
                     its2   LSU anchor, ITS2 >= 100, but ITS1 < 50 or 5.8S < 140
                     other  end-truncated, ITS1-only, undetected
                     chimeric
                     unreadable  characters outside IUPAC (or empty)
  after extraction:  overlong (extracted full ITS > 3,000 bp), not_extracted (pyitsx gave nothing)
  extract          full_ITS span of the full class (pyitsx extract orients 5'->3')
  derep            exact dereplication of the extracted sequences; every input id kept

Everything dropped is counted in classes.tsv / summary.json; nothing is corrected.
Usage: python -m seqhesion intake INPUT.fasta OUTDIR
"""
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

PYITSX = os.path.join(os.path.dirname(sys.executable), 'pyitsx')
RULE = {'ITS1': 50, '5.8S': 140, 'ITS2': 50}
UNANCHORED_MIN = 150                        # ITS1 (no SSU) or ITS2 (no LSU) must be at least this long
IUPAC = set('ACGTRYSWKMBDHVN')
COMP = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
TIP_HEX = 16
MAX_ITS = 3000        # real full ITS tops out near 2,500 bp (Cantharellus); longer is not one ITS


def clean(seq):
    """mm-to-ref's clean_sequence: uppercase; whitespace and gap characters out; no U->T."""
    return re.sub(r'[\s\-.]', '', seq).upper()


def content_id(seq):
    """The input id of a cleaned sequence: orientation-free content hash."""
    return min(hashlib.sha256(seq.encode()).hexdigest(), hashlib.sha256(seq.translate(COMP)[::-1].encode()).hexdigest())


def tip_id(extracted):
    return hashlib.sha256(extracted.upper().encode()).hexdigest()[:TIP_HEX]


class BadInput(ValueError):
    pass


def checked_input(path):
    """[(input_id, cleaned sequence)], refusing ids that are not the content hash, and duplicates."""
    rows, seen, bad = [], set(), []
    for name, seq in read_fasta(path):
        s = clean(seq)
        if name in seen:
            bad.append(f'{name}: duplicate id')
        elif content_id(s) != name:
            bad.append(f'{name}: not the content hash of its sequence')
        seen.add(name)
        rows.append((name, s))
    if bad:
        raise BadInput(f'{len(bad)} of {len(rows)} input records refused, e.g. ' + '; '.join(bad[:5]))
    return rows


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


def run(inputs, outdir, cpus=6, log=print):
    os.makedirs(outdir, exist_ok=True)
    rows_in = checked_input(inputs)
    log(f'{len(rows_in)} input sequences; every id is its content hash')
    readable = os.path.join(outdir, 'readable.fasta')
    unreadable = []
    with open(readable, 'w') as fh:                       # pyitsx rejects anything outside IUPAC
        for name, seq in rows_in:
            if not seq or set(seq) - IUPAC:
                unreadable.append(name)
            else:
                fh.write(f'>{name}\n{seq}\n')
    delim = os.path.join(outdir, 'delimit.tsv')
    if not os.path.exists(delim):
        subprocess.run([PYITSX, 'delimit', '-i', readable, '-o', delim + '.tmp', '--cpus', str(cpus)], check=True)
        os.replace(delim + '.tmp', delim)
    rows = list(csv.DictReader(open(delim), delimiter='\t'))
    cls = {r['seq_id']: classify(r) for r in rows}
    rows += [{'seq_id': n, 'strand': '', 'chimeric': '', **{k: '' for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')}} for n in unreadable]
    cls.update({n: 'unreadable' for n in unreadable})
    missing = [n for n, _ in rows_in if n not in cls]         # readable, but no delimit row at all
    rows += [{'seq_id': n, 'strand': '', 'chimeric': '', **{k: '' for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')}} for n in missing]
    cls.update({n: 'other' for n in missing})
    assert len(cls) == len(rows_in) == len(rows), 'every input gets exactly one class'

    full_in = os.path.join(outdir, 'full_input.fasta')
    with open(full_in, 'w') as fh:
        for name, seq in rows_in:
            if cls.get(name) == 'full':
                fh.write(f'>{name}\n{seq}\n')
    ext = os.path.join(outdir, 'full.extract.fasta')
    if not os.path.exists(ext):
        subprocess.run([PYITSX, 'extract', '-i', full_in, '-o', ext + '.tmp', '--cpus', str(cpus), '--region', 'full_ITS'], check=True)
        os.replace(ext + '.tmp', ext)

    groups = defaultdict(list)
    extracted = set()
    for name, seq in read_fasta(ext):
        name = name.split('|')[0]
        extracted.add(name)
        if len(seq) > MAX_ITS:                 # e.g. a pasted contig that happens to contain an ITS
            cls[name] = 'overlong'
            continue
        groups[seq.upper()].append(name)
    for name, c in cls.items():
        if c == 'full' and name not in extracted:
            cls[name] = 'not_extracted'
    counts = Counter(cls.values())
    counts['full_unanchored'] = sum(1 for r in rows if cls[r['seq_id']] == 'full' and anchors(r) != 'SSU+LSU')
    counts['undetected_by_delimit'] = sum(1 for r in rows if not any(r[k] not in ('', '-', 'None') for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')))
    counts['reverse_strand'] = sum(1 for r in rows if r['strand'] == '-')
    log(f'classes: {dict(counts)}')
    tips = {}
    for seq in groups:
        t = tip_id(seq)
        if t in tips:
            raise RuntimeError(f'tip id collision at {TIP_HEX} hex digits: {t}; lengthen TIP_HEX')
        tips[t] = seq
    reps = sorted(groups.items(), key=lambda kv: (-len(kv[1]), tip_id(kv[0])))
    with open(os.path.join(outdir, 'full.derep.fasta'), 'w') as fh, open(os.path.join(outdir, 'derep_members.tsv'), 'w') as mh:
        mh.write('unique\tsize\tmembers\n')
        for seq, names in reps:
            t = tip_id(seq)
            fh.write(f'>{t}\n{seq}\n')
            mh.write(f'{t}\t{len(names)}\t{",".join(sorted(names))}\n')
    tip_of = {n: tip_id(seq) for seq, names in groups.items() for n in names}
    with open(os.path.join(outdir, 'classes.tsv'), 'w') as fh:
        fh.write('seq_id\tclass\ttip\tstrand\tanchors\tSSU\tITS1\t5.8S\tITS2\tLSU\n')
        for r in rows:
            n = r['seq_id']
            fh.write('\t'.join([n, cls[n], tip_of.get(n, ''), r['strand'], anchors(r) if cls[n] != 'unreadable' else '']
                               + [str(span(r[k])) for k in ('SSU', 'ITS1', '5.8S', 'ITS2', 'LSU')]) + '\n')
    summary = {'input': len(rows), 'classes': dict(counts), 'extracted': sum(len(v) for v in groups.values()),
               'unique': len(groups), 'singletons': sum(len(v) == 1 for v in groups.values()),
               'largest_identical_group': max((len(v) for v in groups.values()), default=0)}
    json.dump(summary, open(os.path.join(outdir, 'summary.json'), 'w'), indent=1)
    log(json.dumps(summary, indent=1))
    return summary


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs')
    ap.add_argument('outdir')
    ap.add_argument('--cpus', type=int, default=6)
    a = ap.parse_args()
    run(a.inputs, a.outdir, a.cpus)
