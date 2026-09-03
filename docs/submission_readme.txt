LaTeX sources for "What Cycle Closure Can and Cannot See in a Relative Binding Free
Energy Network" (Article) and its Supporting Information.

Two documents:
  paper_jctc.tex     the manuscript, which inputs paper_body.tex
  paper_jctc_si.tex  the Supporting Information, which inputs paper_si.tex

They reference each other's figures, tables and sections through the xr package,
which reads the other document's .aux file. Neither is settled until both have been
compiled after the other, so build them in alternation:

  pdflatex paper_jctc      bibtex paper_jctc      pdflatex paper_jctc
  pdflatex paper_jctc_si   bibtex paper_jctc_si   pdflatex paper_jctc_si
  pdflatex paper_jctc      pdflatex paper_jctc_si
  pdflatex paper_jctc      pdflatex paper_jctc_si

A single pass leaves the cross-document numbers stale rather than undefined, which is
why the last passes are repeated. Graphics are in figures/; the preamble sets
\graphicspath accordingly.

Bibliography: refs.bib, cited with natbib \citep as provided by achemso.
