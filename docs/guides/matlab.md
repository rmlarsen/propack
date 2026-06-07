# MATLAB toolbox guide

> This page describes the original **MATLAB** implementation of PROPACK. The
> verbatim original overview is at [`../legacy/Readme.txt`](../legacy/Readme.txt).
> The MATLAB routines ship with detailed on-line help, readable via
> `help <functionname>` at the MATLAB prompt; see also the
> [routine reference](../routines/index.md).

## Package contents

### Main routines

| File | Description |
|------|-------------|
| `lanbpro.m` | Lanczos bidiagonalization with partial reorthogonalization (BPRO). |
| `lansvd.m`  | Sparse SVD routine built on BPRO. |
| `lanpro.m`  | Lanczos tridiagonalization with partial reorthogonalization (PRO). |
| `laneig.m`  | Sparse symmetric eigenvalue routine built on PRO (similar to the `lanso` subroutine in the LANSO package). |

### Computational routines

| File | Description |
|------|-------------|
| `refinebounds.m` | Refines the error bounds based on the gap structure of the Ritz values. |
| `compute_int.m`  | Determines which Lanczos vectors to purge in a reorthogonalization. |
| `reorth.m`       | Iterated (modified or classical) Gram-Schmidt reorthogonalization. |
| `reorth.mex*`    | MEX file using a Fortran 77 version of `reorth.m`. |
| `tqlb.m`         | Eigenvalues and the top and bottom eigenvector elements of a symmetric tridiagonal matrix. |
| `tqlb.mex*`      | MEX files calling a modified version of the EISPACK routine `TQL1`. |
| `bdsqr.m`        | Singular values and the bottom elements of the left singular vectors of a lower bidiagonal matrix. |
| `bdsqr.mex*`     | MEX files calling the LAPACK routine `DBDSQR`. |

### Test scripts and problems

`test.m` (main sparse-SVD test), `testtqlb.m`, and the operator helpers
`Afunc.m` (`Y <- A*X`), `Atfunc.m` (`Y <- A'*X`), `AtAfunc.m` (`Y <- A'*A*X`), and
`Cfunc.m` (`Y <- [0 A; A' 0] * X`). `helio.mat` is a test problem from
helioseismology: a full 212-by-100 matrix with condition number 2.9·10^12. The
`Harwell-Boeing/` directory holds test matrices in MatrixMarket format plus MATLAB
functions for querying, reading, and writing that format.

The MEX versions of `TQLB`, `REORTH`, and `DBDSQR` are included for speed
(compiled for SGI64, SUN, and Linux x86 in the original distribution). Pure-MATLAB
versions are provided for other platforms but can be considerably slower.

## Implementation note: modifications to TQLB

`TQLB` computes the eigenvalues and the top and bottom elements of the eigenvectors
of the Lanczos tridiagonal matrix, which enter into the error bounds. It is based
on an optimized version of the EISPACK routine `TQL1`. Beyond computing the bottom
elements of the normalized eigenvectors, the inner-loop calls to `PYTHAG` were
replaced with in-line code following

> A. Cline and J. Meyering, "Converting EISPACK to run efficiently on a vector
> processor," Tech. Memo., Pleasant Valley Software, Austin TX, 1989.

This decreased execution time by a factor of ~3 versus the original and also
computes the small eigenvalues more accurately. As a test, the eigenvalues of the
standard tridiagonal matrix

```
    (  2  -1            )
    ( -1   2  -1        )
T = (      -1   2   .   )      with  lambda_i = 4*cos(pi/2 * i/(n+1))^2,
    (           .   .  -1)            i = 1, 2, ..., n
    (              -1   2)
```

were computed on an SGI Origin 200 (180 MHz MIPS R10000, 32 KB primary / 1 MB
secondary cache). Timings `T` (seconds) and maximum relative error `E`:

| n | Original TQLB (T) | Original (E) | Modified TQLB (T) | Modified (E) | MATLAB `eig` (T) | `eig` (E) |
|---|------|------|------|------|------|------|
| 100   | 0.020 | 1.4e-12 | 0.007 | 3.1e-13 | 0.024 | 6.2e-14 |
| 1000  | 1.68  | 1.2e-9  | 0.55  | 2.9e-10 | 1.98  | 1.5e-10 |
| 10000 | 154.1 | 1.9e-6  | 52.2  | 8.2e-8  | 181.2 | 3.0e-8  |
