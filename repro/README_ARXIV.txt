ARXIV SUBMISSION PACKAGE
========================

Title: Cutoff and a Gaussian-shift profile for the Rudvalis shuffle
Author: Yunjiang Jiang
Date: October 3, 2026

Suggested arXiv category:
  Primary: math.PR
  Secondary: math.CO

2020 MSC:
  Primary 60J10; secondary 60C05, 20C30.

Suggested comments:
  22 pages. Gives a matching sharp upper bound, total-variation cutoff,
  the n^3 cutoff window, and exact total-variation and chi-square profiles
  for the full labelled Rudvalis shuffle. Includes reproducibility checks.

Main source:
  ../paper.tex

Ancillary reproducibility file:
  checks/verify_rudvalis_profile.py

Build:
  pdflatex ../paper.tex
  pdflatex ../paper.tex

No BibTeX, external figures, custom classes, or nonstandard local style files are required.
The bibliography is included directly in the TeX source.

The finite-dimensional checks are supplementary diagnostics; the manuscript's proof is analytic
and does not depend on them as computer-assisted certificates.
