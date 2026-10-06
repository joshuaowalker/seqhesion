# seqhesion

*Sequence cohesion*: which sequences group together, at every scale the data supports, and how
strongly. Built for fungal ITS barcodes, from many small overlapping trees rather than one large
one.

seqhesion takes a set of sequences and produces a hierarchy of groups, from near-identical
sequences up to relationships across the whole corpus, with the evidence behind every group. It
does not name anything. Labels may travel alongside the input for display, but nothing in the
processing path reads them, and agreement with labels is never a target. Deciding which groups
become named units is curation, and happens downstream.

The result is a map for navigating a large barcode collection: which sequences belong together,
how confidently, and what else is nearby. It is not a phylogeny. Groups at the coarse and corpus
levels are leads for further research, not taxonomic conclusions; ITS alone cannot settle deep
relationships.

## Principles

- **Label-blind.** Names never enter a decision. They may be shown next to the output, and used
  afterwards as a sanity check, never as a target.
- **Many small trees.** No single alignment of a diverse corpus is trustworthy. Each tree covers
  one neighbourhood of 150–180 sequences that align well, and a sequence's place in the hierarchy
  is decided by every tree that holds it.
- **Tree distance and identity distance are not comparable.** Every grouping decision is made
  from the shard trees (model-corrected distance on trimmed columns, sub-resolution branches
  collapsed). Alignment identity (vsearch) is used only to decide what is sampled together: which
  sequences form a component or a shard. The `identity_*` columns of a release are reporting for
  readers, not evidence.
- **Annotate, don't threshold.** Every group carries its confidence measures. Choosing a
  partition is left to the reader.
- **Nothing is corrected, nothing silently dropped.** Sequences are used as given; every input is
  accounted for in a release, with the reason when it was not built.

## How it works

1. **Input.** Every input sequence is identified by a content hash: the smaller of sha256 of the
   cleaned sequence and of its reverse complement. `seqhesion input` prepares any FASTA this way;
   exact copies in either orientation become one input.
2. **Intake.** [pyitsx](https://pypi.org/project/pyitsx/) finds the ITS regions and orients each
   sequence. Sequences with a complete ITS1–5.8S–ITS2 become *tips*, exactly dereplicated (a
   tip's id is a hash of its extracted sequence). ITS2-only, chimeric and other sequences are
   counted by class and not built.
3. **Regions.** Connected components of the graph of 90% centroids, linked at 86% identity. A
   component of 150 or more tips is built on its own. A smaller one draws its rooting context from
   its nearest sequences anywhere in the corpus; that context is pruned before grouping. A lone
   tip gets no tree.
4. **Covers.** Two independent covers of each component by *shards*: a seed and its 149 nearest
   neighbours, topped up so every shard holds at least 30 sequences below 97% identity to its
   seed (150–180 sequences), plus a shared scaffold of up to 30 sequences for rooting. Every
   sequence is in at least three shards of each cover.
5. **Shard trees.** MAFFT L-INS-i, trimAl `-gappyout`, FastTree `-gtr -gamma`; cached under a
   key of exactly what made them. Each tree is rooted on the scaffold, the scaffold pruned, and
   branches shorter than half an expected change collapsed.
6. **Co-association.** For each pair of sequences, every shard holding both records the level at
   which it joins them: the diameter of the smallest clade holding both. The median over shards
   is the pair's join level. Average linkage over observed pairs only gives the hierarchy. Pairs
   no shard holds are not evidence, and parts nothing connects stay a forest. Both covers are
   pooled; a shard drawn by both (same members, so the same tree) counts once.
7. **Evidence.** The hierarchy is cut at levels 0.005 … 0.1 (substitutions per site). Each group
   gets:
   - **cohesion**: the share of shard votes on its pairs that join them at or below the level;
   - **pull**: the share of its members' votes that join them to their strongest outside target;
   - **margin** = cohesion − pull: the confidence measure. It predicts which groups recur when
     the shards are resampled.

   Each group also carries its stem length, spread, nearest outside group, and every shard's
   verdict on it (clade, unresolved or conflict).
8. **Coarse layer** (levels 0.125 … 0.3). Local shards rarely hold two distant sequences, so the
   groups at 0.1 (and lone tips) become *units*. Stratified sparse shards then sample them: one
   random member from each of 150 random units, plus an outgroup from outside the component,
   giving about 24 votes per unit pair. Average linkage over unit-pair medians builds levels
   above the fine ones.
9. **Corpus levels** (0.4, 0.5, 0.75). Each component's groups at 0.3 become super-units, and a
   cover over them, chosen by representative identity, relates components across the corpus. The
   result is one tree over every tip; corpus groups may span many components.
10. **Releases and antenomina.** A release is a set of TSV tables and a JSON manifest. Each group
    has an *antenomen*: a stable id carried from release to release when a group and its
    predecessor are each other's best match and share more than half of the inputs both releases
    hold. `lineage.tsv` logs every match, and an id is never reused.

The fine levels are where the shard trees resolve structure directly (to about 0.075). The coarse
and corpus levels are relative: their groups are stable under resampling, but each is a starting
point for investigation.

## Install

    pip install seqhesion          # or, from a clone: pip install -e '.[test]'

External tools, found on `PATH`: [vsearch](https://github.com/torognes/vsearch),
[MAFFT](https://mafft.cbrc.jp/alignment/software/), [trimAl](https://github.com/inab/trimal) and
[FastTree](http://www.microbesonline.org/fasttree/). To use a particular build, name it in
`SEQHESION_MAFFT`, `SEQHESION_TRIMAL` or `SEQHESION_FASTTREE`. The shard-tree cache is keyed
by each tool's version, so a different version builds new trees rather than reusing old ones.

Alignment is nearly all of the compute. An optimised MAFFT 7.526 whose output is byte-identical to
stock, about 5x faster on Apple silicon and on AVX-512 x86, is at
[joshuaowalker/mafft](https://github.com/joshuaowalker/mafft) (branch `exact-speedups`, with a
technical report in `paper/`). Check any build against stock on a sample of your own shards
before relying on it.

## Use

From a FASTA to a release (every command takes `--procs N`):

    seqhesion input   sequences.fasta INPUT           # content-hash ids, names.tsv, manifest
    seqhesion intake  INPUT/sequences.fasta INTAKE    # pyitsx: extract, orient, dereplicate
    seqhesion regions INTAKE REGIONS                  # writes REGIONS/large.txt and small.txt

    # the fine hierarchy of every component (independent: run them in parallel)
    for r in $(cat REGIONS/large.txt REGIONS/small.txt); do
        seqhesion build $r --out $r/hierarchy.json
    done

    # the coarse layer of every component
    seqhesion sparse-plan INTAKE $(cat REGIONS/large.txt REGIONS/small.txt)
    for r in $(cat REGIONS/large.txt REGIONS/small.txt); do seqhesion sparse-build $r; done

    # corpus levels (corpus-trees --chunk I --of N splits the trees over machines)
    seqhesion corpus-plan   INTAKE CORPUS $(cat REGIONS/large.txt REGIONS/small.txt)
    seqhesion corpus-trees  CORPUS
    seqhesion corpus-layer  CORPUS

    seqhesion release OUT --input INPUT --intake INTAKE --layer --corpus CORPUS \
        $(for r in $(cat REGIONS/large.txt REGIONS/small.txt); do echo $r=$r/hierarchy.json; done) \
        [--previous PREVIOUS_RELEASE]

`build` is `cover`, `prepare` and `hierarchy` in one; they can also be run one at a time. Without
`--corpus`, the coarse layer runs to 0.5 within each component; without `--layer`, a release
stops at 0.1. `--previous` carries group ids forward from an earlier release.

Almost all of the compute is shard trees, one MAFFT L-INS-i alignment and one FastTree per shard.
Trees are cached, so a rebuild after the corpus grows rebuilds only the shards that changed.

## Releases

A release directory holds:

| file | contents |
| --- | --- |
| `manifest.json` | schema (`seqhesion-release/0.5`), seqhesion commit, input fingerprint, method, levels, counts, and the columns of every table |
| `inputs.tsv` | every input: its tip, and whether it was built, not built, or dropped (with its intake class) |
| `tips.tsv`, `tips.fasta` | every tip of a built component, with its extracted, oriented full-ITS sequence |
| `groups.tsv` | one row per distinct group: size, level range, parent, cohesion, pull, margin, stem, spread, nearest group, shard verdicts |
| `group_levels.tsv`, `group_members.tsv` | a group's evidence at each level; its tips |
| `membership.tsv` | per tip and level: its group, how firmly it belongs, its strongest outside pull |
| `lineage.tsv`, `minted.tsv` | how each group id relates to the previous release; every id ever minted |
| `dendrogram/` | the linkage the levels are cut from, as Newick (`corpus.nwk`: one tree of every tip) |
| `shards/` | the prepared shard trees, for drawing and on-demand distances |

The columns are documented in `seqhesion/export.py` and listed in every manifest.

## Status

Version 0.1. seqhesion is in production on a corpus of fungal ITS sequences drawn mainly from
MycoMap, with other sources: about 226,000 input sequences, 147,000 distinct full-ITS tips in
4,172 components. The method was developed and validated in a separate research prototype
(simulation, reproducibility between independent covers, resampling of shards). The draft white
paper in [`docs/whitepaper`](docs/whitepaper) describes the method, its validation and its
limits, and is revised alongside the code.

Not yet in the release: placement of ITS2-only sequences on the hierarchy. The code is in
`seqhesion/place.py`, but running it at corpus scale awaits a batch interface.

## License

BSD 3-Clause; see `LICENSE`.
