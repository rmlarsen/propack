---
title: "PROPACK Fortran 77 guide (v2.1)"
description: "Overview and installation for the classic Fortran 77 PROPACK distribution, version 2.1."
---

# Fortran package guide (classic v2.1)

> This page documents the **classic Fortran 77 distribution (version 2.1, April
> 2005)** as described on the original homepage. For the **current** build (CMake,
> the C++17 API, and the four precision libraries) see the top-level
> [README](https://github.com/rmlarsen/propack/blob/main/README.md). The verbatim original is at
> [`../legacy/README_v2.1.txt`](../legacy/README_v2.1.txt).

## Overview

PROPACK efficiently computes the singular values and singular vectors of a large,
sparse, and/or structured matrix. The basic Krylov-subspace algorithm is Lanczos
bidiagonalization, implemented with partial reorthogonalization. Partial
reorthogonalization often improves performance significantly compared to the
classic Lanczos algorithm with full reorthogonalization; the exact amount depends
on the distribution of the singular values.

Two sets of SVD routines are available, one with and one without implicit
restarting. Implicit restarting computes a given number of singular values and
vectors in a *fixed* amount of memory. The ordinary version uses memory
proportional to the number of iterations required for convergence (generally not
known in advance), but since its total number of matrix-vector multiplications is
usually lower, it can still be the method of choice in many cases.

The main driver routines `DLANSVD` and `DLANSVD_IRL` are found in `dlansvd.F` and
`dlansvd_irl.F`, which also document the input parameters. Example programs for
computing the SVD of sparse matrices in several formats — including the
Harwell-Boeing format — are in each precision's `Examples` directory.

## Installation (classic build)

1. Unpack:

   ```sh
   gunzip PROPACK77.tar.gz
   tar xvf PROPACK77.tar
   ```

2. Edit the make-option file `make.<plat>` in `PROPACK/Make`, where `<plat>` is
   one of `linux_gcc_ia32`, `linux_icc_ia32`, `linux_gcc_ia64`, `linux_icc_ia64`,
   `irix`, `sunos`, or `ibm`. In particular set `LINKFLAGS`, `LINKPATH`, and
   `BLAS` so the BLAS library on your system links correctly. Then run:

   ```sh
   ./configure
   ```

   The `configure` script detects your platform and generates `make.inc`. On
   Intel platforms the optional `-icc` argument selects the Intel `icc`/`ifort`
   toolchain (usually faster than gcc, especially on ia64). On AIX, ia32, ia64,
   and IRIX, passing `-openmp` builds a multi-threaded (OpenMP) version; set
   `OMP_NUM_THREADS` to choose the thread count. The parallelization is very fine
   grained and is mostly suited to large matrices (m, n > 100,000) or shared-memory
   machines with low memory latency.

3. Build:

   ```sh
   make
   ```

   This builds `lib<precision>propack_<PLAT>.a` (the PROPACK routines) and
   `lib<precision>lapack_util_<PLAT>.a` (LAPACK 3.0 routines called by PROPACK),
   where `<precision>` is `s`, `d`, `c`, or `z` (single, double, complex,
   double-complex). Link your program with both, plus your system BLAS. The four
   precision libraries are in the `single`, `double`, `complex8`, and `complex16`
   directories.

## Example programs

For each precision, `example.F` and `example_irl.F` (in `Examples`) illustrate the
non-restarted and implicitly restarted algorithms:

```sh
cd <precision>/Examples
make
example.<PLAT>.x     < example.in
example_irl.<PLAT>.x < example_irl.in
```

The programs read a matrix in Harwell-Boeing format and compute a number of
singular values specified in the input file. Test matrices are provided in
`Examples/illc1850.rra` (single/double) and `Examples/mhd1280b.cua`
(complex8/complex16); more are available from
[MatrixMarket](https://math.nist.gov/MatrixMarket/). The programs can also read
diagonal, coordinate, or dense formats (binary or ASCII).

> **Warning.** Matrices stored in binary format are often incompatible between
> machines of different word size (32- vs 64-bit) or endianness (little-endian
> x86/Itanium/Alpha vs big-endian PPC/Power/MIPS).

## Testing the installation

Reference output is provided in `<precision>/Examples/Output`. From the top-level
directory:

```sh
make; make test; make verify
```

builds the example programs for all precisions, runs them on the provided test
matrices, and verifies the results with `Examples/compare.F`. The error bounds in
the test are generous (meant to catch installation errors). For double and
complex*16 the maximal relative error in the singular values should be ~1e-15;
for single and complex*8, ~1e-6.

## Obtaining BLAS

If your system lacks BLAS, the original guide recommended Kazushige Goto's
high-performance BLAS, or the platform-optimized routines from the ATLAS project
(Netlib). Reference Fortran BLAS and documentation are at
[netlib.org/blas](https://www.netlib.org/blas). (The modern build auto-detects
BLAS/LAPACK and prefers OpenBLAS — see the [README](https://github.com/rmlarsen/propack/blob/main/README.md).)
