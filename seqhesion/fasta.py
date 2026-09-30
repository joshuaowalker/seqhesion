"""FASTA in and out. Ids are the first word of the header."""


def read_fasta(path):
    seqs, sid = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                sid = line[1:].split()[0]
                seqs[sid] = []
            elif sid is not None:
                seqs[sid].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def write_fasta(path, seqs, order=None):
    with open(path, 'w') as f:
        for sid in (order if order is not None else seqs):
            f.write(f'>{sid}\n{seqs[sid]}\n')
