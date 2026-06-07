---
title: "PROPACK history and provenance"
description: "Origin of PROPACK, why the original Stanford homepage is preserved here, and corrections to the original page."
---

# History and provenance

## Origin

PROPACK was written by **Rasmus Munk Larsen**. The work began with his 1998
technical report and Ph.D. research at the Department of Computer Science,
University of Århus, Denmark, and the package was subsequently developed and
distributed while he was at the Solar Oscillations Investigation, W. W. Hansen
Experimental Physics Laboratory (HEPL), Stanford University.

For many years the canonical distribution point was the homepage at
`http://sun.stanford.edu/~rmunk/PROPACK/` (last updated May 16, 2012). That page
hosted the Fortran 77 and MATLAB packages, the technical report, two conference
talks, and a set of test matrices from MatrixMarket.

## Why this archive exists

The original Stanford homepage became effectively invisible to web crawlers and
to AI / document-indexing systems — not because of access rules (its `robots.txt`
never disallowed the PROPACK directory) but because the server is HTTP-only (no
HTTPS), is unreachable from many automated clients, and had essentially no inbound
links. As a result, a substantial body of useful material was at risk of being
lost and was absent from public corpora.

This `docs/` tree preserves that material in clean, openly licensed Markdown and
keeps it alongside the actively maintained source code. Verbatim copies of the
original files are kept under [`legacy/`](legacy/index.md); the rest of the docs
are lightly reformatted for readability.

## Corrections to the original page

- **Report number.** The homepage's Documentation section cited the technical
  report as "DAIMI PB-357." The journal of record (DAIMI Report Series) and its
  DOI give **PB-537** (`10.7146/dpb.v27i537.7070`); the report's title page is
  dated **October 1998** (the homepage said September 1998). The corrected
  citation is used throughout this documentation.
- **Contact address.** The legacy `Readme` listed `rmunk@quake.stanford.edu`;
  the current contact is `rmlarsen@gmail.com`.

## Related distributions and mirrors

PROPACK has been mirrored and wrapped by others, including:

- `github.com/optimizers/PROPACK` — a mirror of the original Stanford site.
- `github.com/scipy/PROPACK` — used by SciPy's sparse SVD routines.
- `github.com/JuliaSmoothOptimizers/PROPACK.jl` — a Julia wrapper.

The repository you are reading is the author's own modern Fortran / C / C++
version.

## License

The software is released under the BSD 3-Clause license. The original 2005 BSD
notice is preserved in [`legacy/license.txt`](legacy/license.txt) and in the
repository [`LICENSE`](https://github.com/rmlarsen/propack/blob/main/LICENSE). Documentation is under CC BY 4.0; see
[`LICENSE-docs.md`](LICENSE-docs.md).
