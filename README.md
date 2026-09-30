# seqhesion

*Sequence cohesion*: which sequences group together, at every scale the data supports, and how
strongly. Built for fungal ITS, from many small overlapping trees rather than one large one.

seqhesion takes sequences and produces a hierarchy of groups, each with the evidence behind it.
It does not name anything. Labels may travel alongside the input for display, but nothing in the
processing path reads them, and agreement with labels is never a target. Choosing which groups
become named units is curation, and happens downstream (for MycoMap, in mm-to-ref); curated
labels come back only as test cases.

## How it works

1. **Intake.** Orientation and the full-ITS rule (pyitsx), exact dereplication. Nothing is
   corrected; everything dropped is counted.
2. **Regions.** Connected components of a 90%-centroid similarity graph, built independently.
3. **Covers.** Two independent covers of each region by *shards*: a seed and its nearest
   neighbours (150), topped up so each shard holds at least 30 sequences below 97% identity to
   its seed, plus a shared scaffold for rooting. Every sequence is in at least three shards.
4. **Shard trees.** MAFFT L-INS-i, trimAl -gappyout, FastTree -gtr -gamma; cached, keyed to
   exactly what made them. Each is rooted on the scaffold, the scaffold pruned, and branches
   under half an expected change collapsed.
5. **Co-association.** For each pair of sequences, every shard holding both records the level at
   which it joins them (the diameter of the smallest clade holding both; inside a polytomy,
   sequences the data cannot tell apart join at their own span). The median over shards is the
   pair's join level. Average linkage over *observed pairs only* gives one hierarchy; pairs no
   shard holds are no evidence, and parts nothing connects stay a forest.
6. **Evidence.** The hierarchy is cut at levels 0.005 ... 0.1. Per group: cohesion (shard votes
   joining its pairs), replication (the same votes from the other cover, which the hierarchy
   never saw), held (share of its pairs any shard holds). Per sequence: membership and the
   strongest outside pull. Per distinct group: every shard's verdict (clade / unresolved /
   conflict), stem length, spread, nearest outside group and the gap to it.

The shards support structure to about 0.075; coarser levels are relative.

## Use

    python -m seqhesion intake    INPUT.fasta INTAKE_DIR
    python -m seqhesion regions   INTAKE_DIR REGIONS_DIR
    python -m seqhesion build     REGION_DIR --out hierarchy.json     # cover + prepare + hierarchy

or the steps of `build` one at a time: `cover`, `prepare`, `hierarchy`.

External tools: vsearch, MAFFT, trimAl, FastTree, pyitsx. Python: numpy, scipy, ete3.

## Status

Moved out of the *ubertree* research prototype on 2026-09-30, where the method was developed
and validated (simulation, cross-cover reproducibility, the stitch it replaced). The history and
the experiments stay there. Reproduces that prototype's hierarchies byte for byte.

Next: an export for downstream consumers (TSV tables and a JSON manifest), sequence ids that
are content hashes (intake still numbers unique sequences by abundance), release-to-release
stable group ids, a compute plan for the full corpus.
