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
    own = set(seqs)
    small = json.load(open(rd / 'region.json')).get('kind') == 'small'
    nbrs = cover.load_knn(rd / 'knn.tsv')
    if small:            # neighbours from the whole corpus; seeds only from the component's own tips
        seqs.update(read_fasta(rd / 'context.fasta'))
        nbrs = {t: [x for x in nbrs.get(t, []) if x != t] for t in sorted(own)}
    ident = cover.load_identity(rd / 'knn.tsv') if a.quota else None
    man = cover.plan(nbrs, own, cover.scaffold_of(rd / 'comp_cent_0.85.fasta'), 2, a.depth, a.size,
                     a.quota, a.quota_id, ident, log=_log)
    if small:
        man = cover.as_context(man, own)
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
    res = hierarchy.build(rd / a.covers, tips, levels, a.join, a.procs, log=_log)
    res.pop('shard_trees')
    json.dump(res, open(a.out, 'w'), separators=(',', ':'))
    _log(f'wrote {a.out}')


def cmd_build(a):
    a.covers_dir = str(Path(a.region) / a.covers)
    cmd_cover(a)
    cmd_prepare(a)
    cmd_hierarchy(a)


def cmd_sparse_plan(a):
    """Sparse shards for the large components among REGIONS (small ones need none: sparse-build
    makes their layer from their local trees). Needs the corpus, for the outgroups: run where the
    intake is."""
    from . import sparse
    from .fasta import read_fasta
    large = []
    for r in a.regions:
        rd = Path(r)
        if json.load(open(rd / 'region.json')).get('kind') == 'small':
            continue
        large.append(rd)
    og = sparse.outgroups(large, a.intake, a.procs, log=_log)
    corpus = read_fasta(Path(a.intake) / 'full.derep.fasta')
    for rd in large:
        h = json.load(open(rd / a.hierarchy))
        k = sparse.plan(rd, h, og[rd.name], corpus)
        _log(f'{rd.name}: {k} sparse shards' if k else f'{rd.name}: fewer than 3 units, no sparse shards')


def cmd_sparse_build(a):
    """The coarse layer of one region: from its planned sparse shards (trees built here), else from
    its local trees (small components, and large ones with fewer than 3 units)."""
    from . import sparse
    rd = Path(a.region)
    if (rd / 'sparse' / 'manifest.json').exists():
        sparse.build(rd, a.procs, log=_log)
        return
    h = json.load(open(rd / a.hierarchy))
    _, members = sparse.units(h)
    if len(members) < 2:
        (rd / 'sparse').mkdir(exist_ok=True)
        sparse._write(rd / 'sparse' / 'layer.json', sparse.empty_layer(h))
        _log(f'{rd.name}: one unit, empty layer')
        return
    sparse.build_local(rd, h, a.covers, a.procs, log=_log)


def cmd_corpus_plan(a):
    """Corpus levels, step 1 (needs the corpus): super-units, their neighbours, the cover over them."""
    from . import corpus
    corpus.plan(a.regions, a.intake, a.out, a.hierarchy, a.depth, threads=a.procs, log=_log, seed=a.seed)


def cmd_corpus_trees(a):
    """Corpus levels, step 2: the trees of one chunk of the plan (chunks are independent: AWS array)."""
    from . import corpus
    corpus.trees(a.out, a.chunk, a.of, a.procs, log=_log)


def cmd_corpus_layer(a):
    """Corpus levels, step 3: prepare every tree, build OUT/layer.json."""
    from . import corpus
    corpus.layer(a.out, a.procs, log=_log)


def cmd_release(a):
    from . import export, regions
    if a.out == 'auto':
        if not a.releases:
            raise SystemExit("OUT_DIR 'auto' needs --releases DIR")
        a.out = str(Path(a.releases) / export.next_name(a.releases, 'f'))
        _log(f'release name: {Path(a.out).name}')
    from . import corpus as corpus_mod, sparse
    if a.corpus and not a.layer:
        raise SystemExit('--corpus needs --layer (its units are the coarse groups at 0.3)')
    corpus = None
    if a.corpus:
        cd = Path(a.corpus)
        corpus = (json.load(open(cd / 'units.json')), json.load(open(cd / 'layer.json')))
    comps = []
    for spec in a.components:
        rd, hj = spec.split('=', 1)
        h = json.load(open(hj))
        if a.layer:
            lj = Path(rd) / 'sparse' / 'layer.json'
            if not lj.exists():
                raise SystemExit(f'{rd}: --layer, but no sparse/layer.json (run sparse-plan / sparse-build)')
            h = sparse.add_layer(h, json.load(open(lj)), max_level=corpus_mod.UNIT_LEVEL if a.corpus else None)
        comps.append((rd, h))
    export.write_release(a.out, a.input, a.intake, comps, regions.components(a.intake, a.min_id), a.procs, log=_log,
                         settings={'component_min_id': a.min_id}, previous=a.previous, place_its2=not a.no_place_its2,
                         corpus=corpus)
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
        # always two covers, drawn independently; the hierarchy pools them (cross-cover twins once)
        p.add_argument('--depth', type=int, default=3)
        p.add_argument('--size', type=int, default=150)
        p.add_argument('--quota', type=int, default=30,
                       help='every shard gets at least this many members below --quota-id identity to its seed')
        p.add_argument('--quota-id', type=float, default=97.0)

    def hier_args(p):
        p.add_argument('--levels', default=','.join(str(x) for x in (0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1)))
        p.add_argument('--join', default='identical', choices=('diameter', 'compatible', 'identical'))

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
    p = add('sparse-plan', cmd_sparse_plan, 'intake')
    p.add_argument('regions', nargs='+')
    p.add_argument('--hierarchy', default='hierarchy.json', help="the fine hierarchy's file name in each region dir")
    p = add('sparse-build', cmd_sparse_build, 'region')
    p.add_argument('--hierarchy', default='hierarchy.json')
    p.add_argument('--covers', default='covers')
    p = add('corpus-plan', cmd_corpus_plan, 'intake', 'out')
    p.add_argument('regions', nargs='+')
    p.add_argument('--hierarchy', default='hierarchy.json')
    p.add_argument('--depth', type=int, default=6)
    p.add_argument('--seed', type=int, default=0, help='an independent plan, for a stability check')
    p = add('corpus-trees', cmd_corpus_trees, 'out')
    p.add_argument('--chunk', type=int, default=0)
    p.add_argument('--of', type=int, default=1)
    add('corpus-layer', cmd_corpus_layer, 'out')
    p = add('release', cmd_release, 'out')
    p.add_argument('--releases', help="with OUT_DIR 'auto': the releases directory; the release is named YYYYMMDD.NN + f")
    p.add_argument('--input', required=True, help='the input directory (its manifest.json is recorded)')
    p.add_argument('--intake', required=True)
    p.add_argument('components', nargs='+', help='REGION_DIR=HIERARCHY.json')
    p.add_argument('--min-id', type=float, default=80.0, help='the centroid identity components were cut at')
    p.add_argument('--previous', help='the previous release directory, whose group ids (antenomina) are carried forward')
    p.add_argument('--latest', action='store_true', help="point <releases>/latest at this release")
    p.add_argument('--layer', action='store_true',
                   help="add each region's coarse layer (sparse/layer.json, from sparse-plan / sparse-build) above 0.1")
    p.add_argument('--corpus', help='a corpus layer directory (corpus-plan / -trees / -layer): its levels follow the '
                                    "components', whose layers then stop at 0.3")
    p.add_argument('--no-place-its2', action='store_true',
                   help='skip ITS2-only placement (hours at scale on one machine); those inputs stay dropped:its2')
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == '__main__':
    main()
