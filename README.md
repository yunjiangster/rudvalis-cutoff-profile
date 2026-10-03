# Cutoff and a Gaussian-shift profile for the Rudvalis shuffle

This repository contains the manuscript and reproducibility material for:

**Yunjiang Jiang, _Cutoff and a Gaussian-shift profile for the Rudvalis shuffle_ (2026).**

The Rudvalis shuffle removes the top card and reinserts it uniformly at random in one of the bottom two positions. The manuscript proves a sharp total-variation cutoff for the full labelled shuffle at

[
\frac{n^3\log n}{8\pi^2}
]

with an order-(n^3) cutoff window. More precisely, for fixed (s\in\mathbb R), at

[
t_n(s)=\left\lfloor \frac{n^3}{8\pi^2}(\log n+2s)\right\rfloor,
]

the claimed limiting profiles are

[
d_n(t_n(s))\to \operatorname{erf}(e^{-s}/2),
\qquad
\chi_n^2(t_n(s))\to \exp(2e^{-2s})-1.
]

## Files

- `paper.pdf` — compiled manuscript.
- `paper.tex` — self-contained LaTeX source.
- `repro/` — arXiv/reproducibility source package and verification material from the research writeup.

## Status

The manuscript records a complete analytic derivation developed and audited in-session. It has not yet been independently refereed or formally verified. The finite computations in `repro/` are diagnostic checks of identities, constants, and normalizations; they are not substitutes for the analytic proof.

## Author

Yunjiang Jiang
