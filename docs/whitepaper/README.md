# seqhesion white paper (draft)

`seqhesion.tex` describes the method, its implementation and the lab results behind it.
`references.bib` holds its references; every entry carries a DOI where one exists.

Build (needs a TeX distribution with latexmk, pdflatex, bibtex, natbib and TikZ):

    latexmk -pdf seqhesion.tex

Clean up intermediate files with `latexmk -c` (or `latexmk -C` to remove the PDF too).

The paper is a companion to this repository, revised in place as the method changes. Red
**[TODO: ...]** markers in the PDF are placeholders still to be filled in.
