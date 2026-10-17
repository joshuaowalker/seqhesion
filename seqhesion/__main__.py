"""seqhesion: a label-blind hierarchy of sequence groups from many small overlapping trees.

  seqhesion input        FASTA INPUT_DIR                 content-hash ids, names.tsv, manifest.json
  seqhesion intake       INPUT_DIR/sequences.fasta INTAKE_DIR
  seqhesion regions      INTAKE_DIR REGIONS_DIR          writes REGIONS_DIR/large.txt and small.txt
  seqhesion build        REGION_DIR --out HIERARCHY.json  (= cover + prepare + hierarchy)
  seqhesion cover        REGION_DIR [--reuse OLD_REGIONS_DIR] [--plan-only]
  seqhesion cover-trees-plan REGION_LIST JOBS.json        trees the planned covers still need
  seqhesion cover-trees  JOBS.json REGIONS_ROOT [--chunk I --of N]   build them (e.g. an AWS array)
  seqhesion prepare      REGION_DIR/covers
  seqhesion hierarchy    REGION_DIR --out HIERARCHY.json
  seqhesion stem-nodes   REGION_DIR --hierarchy H.json --out H2.json   (stem nodes for an older hierarchy)
  seqhesion sparse-plan  INTAKE_DIR REGION_DIR ...       coarse layer (levels 0.125 ...)
  seqhesion sparse-build REGION_DIR
  seqhesion corpus-plan  INTAKE_DIR CORPUS_DIR REGION_DIR ...   corpus levels (0.4, 0.5, 0.75)
  seqhesion corpus-trees CORPUS_DIR [--chunk I --of N]
  seqhesion corpus-layer CORPUS_DIR
  seqhesion corpus-verify CORPUS_DIR --round R     then corpus-trees --round R, corpus-layer; repeat
  seqhesion release      OUT_DIR|auto [--releases DIR] --input INPUT_DIR --intake INTAKE_DIR
                         [--layer [--corpus CORPUS_DIR]] [--previous RELEASE] REGION_DIR=HIERARCHY.json ...
"""
import argparse
import json
import sys
from pathlib import Path


def _log(msg):
    print(msg, flush=True)


def cmd_input(a):
    from . import inputs
    inputs.prepare(a.fasta, a.out, log=_log)


def cmd_intake(a):
    from . import intake
    intake.run(a.input, a.outdir, a.procs, log=_log)


def cmd_regions(a):
    from . import regions
    regions.regions(a.intake, a.out, a.min_id, threads=a.procs, log=_log)


def cmd_cover(a):
    from . import cover
    rd = Path(a.region)
    cd = rd / a.covers
    cd.mkdir(exist_ok=True)
    man, seqs = cover.plan_region(rd, a.depth, a.size, a.quota, a.quota_id, log=_log)
    json.dump(man, open(cd / 'manifest.json', 'w'))
    reuse = [str(Path(r) / rd.name / a.covers / 'trees') for r in a.reuse]
    if getattr(a, 'plan_only', False):
        if reuse:
            cover.link_reusable(man, seqs, cd / 'trees', reuse, log=_log)
        return
    cover.build_trees(man, seqs, cd / 'trees', a.procs, log=_log, reuse=reuse)


def cmd_cover_trees_plan(a):
    """The trees the planned covers of the listed regions still need, as a job list (JSON) for
    cover-trees: computed where the earlier trees are, so a remote chunk needs none of them."""
    from . import cover
    regs = [l.strip() for l in open(a.regions) if l.strip()]
    jobs = cover.missing_trees(regs, a.covers, log=_log)
    json.dump({'covers': a.covers, 'jobs': jobs}, open(a.out, 'w'))


def cmd_cover_trees(a):
    """One chunk of a cover-trees-plan job list (chunks are independent: an AWS array). Needs each
    listed region's comp.fasta (+ context.fasta) under REGIONS_ROOT; writes its trees there."""
    from . import cover
    j = json.load(open(a.jobs))
    cover.build_jobs(j['jobs'], a.regions_root, j['covers'], a.chunk, a.of, a.procs, log=_log)


def cmd_place_plan(a):
    """ITS2-only placement, step 1 (where the intake is): statuses that need no tree, and one job per
    (query, shard) for the rest, against the components REGIONS will release (same --min-id)."""
    from . import export, place, regions
    from .fasta import read_fasta
    classes = export.read_tsv(Path(a.intake) / 'classes.tsv')
    comp_of = export.component_of(classes, regions.components(a.intake, a.min_id))
    covers_of = {}
    for r in a.regions:
        rd = Path(r)
        covers_of[comp_of[next(iter(read_fasta(rd / 'comp.fasta')))]] = rd / a.covers
    plan = place.plan_inputs(a.intake, classes, comp_of, covers_of, a.procs, log=_log)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    json.dump(plan, open(out / 'plan.json.tmp', 'w'), separators=(',', ':'))
    (out / 'plan.json.tmp').replace(out / 'plan.json')


def cmd_place_shards(a):
    """ITS2-only placement, step 2: one chunk of the plan's shard re-inferences (an AWS array)."""
    from . import place
    place.place_shards(a.out, a.regions_root, a.chunk, a.of, a.procs, log=_log)


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


def cmd_stem_nodes(a):
    """Add stem nodes to a hierarchy built before they existed (no tree is built)."""
    from . import hierarchy
    rd = Path(a.region)
    h = json.load(open(rd / a.hierarchy))
    hierarchy.add_stem_nodes(h, rd / a.covers, a.procs, log=_log)
    tmp = Path(f'{a.out}.tmp')
    json.dump(h, open(tmp, 'w'), separators=(',', ':'))
    tmp.replace(a.out)
    _log(f'wrote {a.out}')


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
        sparse.build(rd, a.procs, log=_log, reuse=[str(Path(r) / rd.name / 'sparse' / 'trees') for r in a.reuse])
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
    """Corpus levels, step 2: the trees of one chunk of the plan (chunks are independent: AWS array);
    with --round, only that verification round's shards."""
    from . import corpus
    corpus.trees(a.out, a.chunk, a.of, a.procs, log=_log, round_=a.round, reuse=a.reuse)


def cmd_corpus_verify(a):
    """Corpus levels, verification round R: shards for every merge resting on too few votes, appended
    to the plan; then corpus-trees --round R and corpus-layer. Exit status 3 when nothing is thin."""
    from . import corpus
    if corpus.verify(a.out, a.round, a.procs, log=_log, intake_dir=a.intake) == 0:
        raise SystemExit(3)


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
    if a.placements and a.no_place_its2:
        raise SystemExit('--placements and --no-place-its2 contradict each other')
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
                         corpus=corpus, placements=a.placements)
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

    add('input', cmd_input, 'fasta', 'out')
    add('intake', cmd_intake, 'input', 'outdir')
    p = add('regions', cmd_regions, 'intake', 'out')
    p.add_argument('--min-id', type=float, default=86.0, help='identity (%%) at which two 90%% centroids link components')

    def cover_args(p):
        # always two covers, drawn independently; the hierarchy pools them (cross-cover twins once)
        p.add_argument('--reuse', action='append', default=[], metavar='REGIONS_DIR',
                       help="an earlier cycle's regions dir: any of its trees whose key is exactly a shard's is "
                            "hard-linked (by content, whatever its shard name); repeatable")
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
    p.add_argument('--plan-only', action='store_true', help='write the plan (and link reusable trees), build nothing')
    cover_args(p)
    p = add('cover-trees-plan', cmd_cover_trees_plan, 'regions')
    p.add_argument('out')
    p.add_argument('--covers', default='covers')
    p = add('cover-trees', cmd_cover_trees, 'jobs')
    p.add_argument('regions_root')
    p.add_argument('--chunk', type=int, default=0)
    p.add_argument('--of', type=int, default=1)
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
    p = add('stem-nodes', cmd_stem_nodes, 'region')
    p.add_argument('--hierarchy', default='hierarchy.json')
    p.add_argument('--covers', default='covers')
    p.add_argument('--out', required=True)
    p = add('sparse-plan', cmd_sparse_plan, 'intake')
    p.add_argument('regions', nargs='+')
    p.add_argument('--hierarchy', default='hierarchy.json', help="the fine hierarchy's file name in each region dir")
    p = add('sparse-build', cmd_sparse_build, 'region')
    p.add_argument('--hierarchy', default='hierarchy.json')
    p.add_argument('--covers', default='covers')
    p.add_argument('--reuse', action='append', default=[], metavar='REGIONS_DIR',
                   help="an earlier cycle's regions dir: sparse trees whose key is exactly a shard's are hard-linked")
    p = add('corpus-plan', cmd_corpus_plan, 'intake', 'out')
    p.add_argument('regions', nargs='+')
    p.add_argument('--hierarchy', default='hierarchy.json')
    p.add_argument('--depth', type=int, default=6)
    p.add_argument('--seed', type=int, default=0, help='an independent plan, for a stability check')
    p = add('corpus-trees', cmd_corpus_trees, 'out')
    p.add_argument('--chunk', type=int, default=0)
    p.add_argument('--of', type=int, default=1)
    p.add_argument('--round', type=int, help='only the shards of this verification round')
    p.add_argument('--reuse', action='append', default=[], metavar='TREES_DIR',
                   help="an earlier corpus plan's trees dir: trees whose key is exactly a shard's are hard-linked")
    p = add('corpus-verify', cmd_corpus_verify, 'out')
    p.add_argument('--round', type=int, required=True)
    p.add_argument('--intake', help="the intake dir, when the plan does not record it (plans before 2026-10-16)")
    add('corpus-layer', cmd_corpus_layer, 'out')
    p = add('place-plan', cmd_place_plan, 'intake', 'out')
    p.add_argument('regions', nargs='+', help="the region dirs the release will build (the release's specs)")
    p.add_argument('--min-id', type=float, default=86.0, help="as the release's --min-id")
    p.add_argument('--covers', default='covers')
    p = add('place-shards', cmd_place_shards, 'out')
    p.add_argument('regions_root', help='the directory holding the planned region dirs (by name)')
    p.add_argument('--chunk', type=int, default=0)
    p.add_argument('--of', type=int, default=1)
    p = add('release', cmd_release, 'out')
    p.add_argument('--releases', help="with OUT_DIR 'auto': the releases directory; the release is named YYYYMMDD.NN + f")
    p.add_argument('--input', required=True, help='the input directory (its manifest.json is recorded)')
    p.add_argument('--intake', required=True)
    p.add_argument('components', nargs='+', help='REGION_DIR=HIERARCHY.json')
    p.add_argument('--min-id', type=float, default=86.0, help='the centroid identity components were cut at (as for regions)')
    p.add_argument('--previous', help='the previous release directory, whose group ids (antenomina) are carried forward')
    p.add_argument('--latest', action='store_true', help="point <releases>/latest at this release")
    p.add_argument('--layer', action='store_true',
                   help="add each region's coarse layer (sparse/layer.json, from sparse-plan / sparse-build) above 0.1")
    p.add_argument('--corpus', help='a corpus layer directory (corpus-plan / -trees / -layer): its levels follow the '
                                    "components', whose layers then stop at 0.3")
    p.add_argument('--no-place-its2', action='store_true',
                   help='skip ITS2-only placement (hours at scale on one machine); those inputs stay dropped:its2')
    p.add_argument('--placements', metavar='PLACE_DIR',
                   help='finish a placement run (place-plan, then every place-shards chunk) instead of placing here')
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == '__main__':
    main()
