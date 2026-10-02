"""seqhesion: a label-blind hierarchy of sequence groups from many small overlapping trees.

  python -m seqhesion intake    INPUT.fasta INTAKE_DIR
  python -m seqhesion regions   INTAKE_DIR REGIONS_DIR
  python -m seqhesion cover     REGION_DIR [--covers covers] [--quota 30]
  python -m seqhesion prepare   COVERS_DIR
  python -m seqhesion hierarchy REGION_DIR [--covers covers] --out hierarchy.json
  python -m seqhesion build     REGION_DIR [--covers covers] --out hierarchy.json   (cover + prepare + hierarchy)
  python -m seqhesion release   OUT_DIR|auto [--releases DIR] --input INPUT_DIR --intake INTAKE_DIR REGION_DIR=HIERARCHY.json ...
"""
import argparse
import json
import sys
from pathlib import Path


def _log(msg):
    print(msg, flush=True)


def cmd_intake(a):
    from . import intake
    intake.run(a.input, a.outdir, a.procs, log=_log)


def cmd_regions(a):
    from . import regions
    regions.regions(a.intake, a.out, a.min_id, threads=a.procs, log=_log)


def cmd_cover(a):
    from . import cover
    from .fasta import read_fasta
    rd = Path(a.region)
    cd = rd / a.covers
    cd.mkdir(exist_ok=True)
    seqs = read_fasta(rd / 'comp.fasta')
    nbrs = cover.load_knn(rd / 'knn.tsv')
    ident = cover.load_identity(rd / 'knn.tsv') if a.quota else None
    man = cover.plan(nbrs, set(seqs), cover.scaffold_of(rd / 'comp_cent_0.85.fasta'), 2, a.depth, a.size,
                     a.quota, a.quota_id, ident, log=_log)
    json.dump(man, open(cd / 'manifest.json', 'w'))
    cover.build_trees(man, seqs, cd / 'trees', a.procs, log=_log)


def cmd_prepare(a):
    from . import shardtrees
    for c in (0, 1):
        shardtrees.export(a.covers_dir, c, a.procs, log=_log)


def cmd_hierarchy(a):
    from . import hierarchy
    from .fasta import read_fasta
    rd = Path(a.region)
    tips = sorted(read_fasta(rd / 'comp.fasta'))
    levels = [float(x) for x in a.levels.split(',')]
    res = hierarchy.build(rd / a.covers, tips, levels, a.join, a.primary, a.procs, log=_log)
    res.pop('shard_trees')
    json.dump(res, open(a.out, 'w'), separators=(',', ':'))
    _log(f'wrote {a.out}')


def cmd_build(a):
    a.covers_dir = str(Path(a.region) / a.covers)
    cmd_cover(a)
    cmd_prepare(a)
    cmd_hierarchy(a)


def cmd_release(a):
    from . import export, regions
    if a.out == 'auto':
        if not a.releases:
            raise SystemExit("OUT_DIR 'auto' needs --releases DIR")
        a.out = str(Path(a.releases) / export.next_name(a.releases, 'f'))
        _log(f'release name: {Path(a.out).name}')
    comps = []
    for spec in a.components:
        rd, hj = spec.split('=', 1)
        comps.append((rd, json.load(open(hj))))
    export.write_release(a.out, a.input, a.intake, comps, regions.components(a.intake, a.min_id), a.procs, log=_log,
                         settings={'component_min_id': a.min_id}, previous=a.previous)
    if a.latest:
        link = Path(a.out).parent / 'latest'
        tmp = Path(a.out).parent / '.latest.tmp'
        if tmp.is_symlink():
            tmp.unlink()
        tmp.symlink_to(Path(a.out).name)
        tmp.replace(link)
        _log(f'latest -> {Path(a.out).name}')


def main(argv=None):
    ap = argparse.ArgumentParser(prog='seqhesion', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def add(name, fn, *pos):
        p = sub.add_parser(name)
        for x in pos:
            p.add_argument(x)
        p.add_argument('--procs', type=int, default=8)
        p.set_defaults(fn=fn)
        return p

    add('intake', cmd_intake, 'input', 'outdir')
    p = add('regions', cmd_regions, 'intake', 'out')  # --min-id: identity (%) linking 90% centroids
    p.add_argument('--min-id', type=float, default=80.0)

    def cover_args(p):
        # always two covers: the hierarchy is built from one, and the other replicates it
        p.add_argument('--depth', type=int, default=3)
        p.add_argument('--size', type=int, default=150)
        p.add_argument('--quota', type=int, default=30,
                       help='every shard gets at least this many members below --quota-id identity to its seed')
        p.add_argument('--quota-id', type=float, default=97.0)

    def hier_args(p):
        p.add_argument('--levels', default=','.join(str(x) for x in (0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1)))
        p.add_argument('--join', default='identical', choices=('diameter', 'compatible', 'identical'))
        p.add_argument('--primary', type=int, default=0, choices=(0, 1))

    p = add('cover', cmd_cover, 'region')
    p.add_argument('--covers', default='covers')
    cover_args(p)
    add('prepare', cmd_prepare, 'covers_dir')
    p = add('hierarchy', cmd_hierarchy, 'region')
    p.add_argument('--covers', default='covers')
    p.add_argument('--out', required=True)
    hier_args(p)
    p = add('build', cmd_build, 'region')
    p.add_argument('--covers', default='covers')
    p.add_argument('--out', required=True)
    cover_args(p)
    hier_args(p)
    p = add('release', cmd_release, 'out')
    p.add_argument('--releases', help="with OUT_DIR 'auto': the releases directory; the release is named YYYYMMDD.NN + f")
    p.add_argument('--input', required=True, help='the input directory (its manifest.json is recorded)')
    p.add_argument('--intake', required=True)
    p.add_argument('components', nargs='+', help='REGION_DIR=HIERARCHY.json')
    p.add_argument('--min-id', type=float, default=80.0, help='the centroid identity components were cut at')
    p.add_argument('--previous', help='the previous release directory, whose group ids (antenomina) are carried forward')
    p.add_argument('--latest', action='store_true', help="point <releases>/latest at this release")
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == '__main__':
    main()
