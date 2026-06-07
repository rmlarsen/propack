---
title: "PROPACK technical report and talks"
description: "The 1998 technical report (DAIMI PB-537) introducing BPRO, with its abstract, plus two conference talks."
---

# Technical report and talks

## Technical report

### Lanczos bidiagonalization with partial reorthogonalization

**Rasmus Munk Larsen.** Department of Computer Science, University of Århus,
Ny Munkegade building 540, DK-8000 Århus C, Denmark. October 1998.
Technical Report **DAIMI PB-537**.
DOI: [10.7146/dpb.v27i537.7070](https://doi.org/10.7146/dpb.v27i537.7070)

Local copies: [PDF](paper.pdf) · [gzipped PostScript](paper.ps.gz).
This report describes the work that led to PROPACK.

#### Abstract

A partial reorthogonalization procedure (BPRO) for maintaining semi-orthogonality
among the left and right Lanczos vectors in the Lanczos bidiagonalization (LBD)
is presented. The resulting algorithm is mathematically equivalent to the
symmetric Lanczos algorithm with partial reorthogonalization (PRO) developed by
Simon, but works directly on the Lanczos bidiagonalization of *A*. For computing
the singular values and vectors of a large sparse matrix with high accuracy, the
BPRO algorithm uses only half the amount of storage and a factor of 3–4 less
work compared to methods based on PRO applied to an equivalent symmetric system.
Like PRO, the algorithm presented here is based on simple recurrences, which
enable it to monitor the loss of orthogonality among the Lanczos vectors directly
without forming inner products. These recurrences are used to develop a Lanczos
bidiagonalization algorithm with partial reorthogonalization, which has been
implemented in a MATLAB package for sparse SVD and eigenvalue problems called
PROPACK. Numerical experiments with the routines from PROPACK are conducted using
a test problem from inverse helioseismology to illustrate the properties of the
method. In addition, a number of test matrices from the Harwell-Boeing collection
are used to compare the accuracy and efficiency of the MATLAB implementations of
BPRO and PRO with the `svds` routine in MATLAB 5.1, which uses an implicitly
restarted Lanczos algorithm.

#### Contents

1. Introduction (notation)
2. Lanczos bidiagonalization
3. The Lanczos algorithm and sparse SVD calculations — equivalent symmetric
   problems; fundamental error analysis for SVD calculations; the Lanczos
   algorithm; approximate SVD using Lanczos on the matrix *C*; connection to
   Lanczos bidiagonalization
4. Sparse least squares — least squares solvers based on LBD; algorithms for
   problems with multiple right-hand sides
5. Lanczos bidiagonalization in finite precision arithmetic — Lanczos with no
   reorthogonalization; stabilizing Lanczos using reorthogonalization
6. A partial reorthogonalization algorithm for LBD — tracking the loss of
   orthogonality; computing the reorthogonalization; computing small singular
   values; a correction to the LANSO package; a hybrid method to improve
   performance on cache-based architectures
7. Numerical experiments — accuracy of the computed singular values; efficiency
   of the reorthogonalization methods; comparison of sparse SVD algorithms
8. Conclusion
   - References
   - Appendix A: Singular values of the test matrices
   - Appendix B: PROPACK — sparse SVD and eigenvalue routines in MATLAB

> **Note on the report number.** The original homepage listed this report as
> "DAIMI PB-357, September 1998." The journal of record (DAIMI Report Series) and
> its DOI give **PB-537**, and the report's own title page is dated **October
> 1998**; those are used here.

## Talks

### Computing the SVD for large and sparse matrices

R. M. Larsen, SCCM, Stanford University. Presented June 16, 2000 at Lawrence
Berkeley National Laboratory. [PDF](talk-lbnl-2000.pdf).

### Combining implicit restart and partial reorthogonalization in Lanczos bidiagonalization

R. M. Larsen, SCCM, Stanford University. Presented April 18, 2001 at UC Berkeley.
A presentation of the implicitly restarted SVD algorithm available in PROPACK
version 1.0, including performance measurements of the Fortran version and
comparisons with PLANSO and ARPACK. [PDF](talk-berkeley-2001.pdf).
