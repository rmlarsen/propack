# PROPACK

**Software for computing the singular value decomposition (SVD) of large and
sparse or structured matrices**, based on Lanczos bidiagonalization with partial
reorthogonalization (BPRO).

PROPACK contains a set of routines for computing the SVD of large and sparse or
structured matrices. The SVD routines are based on the Lanczos bidiagonalization
algorithm with partial reorthogonalization (BPRO). The Lanczos routines can also
be used directly, and form the basis of efficient algorithms for solving linear
systems of equations and linear least squares problems, in particular for
systems with multiple right-hand sides.

Partial reorthogonalization often improves performance significantly compared to
the classic Lanczos algorithm with full reorthogonalization; the exact amount of
improvement depends on the distribution of the singular values.

> **History.** PROPACK was originally developed by Rasmus Munk Larsen and
> distributed from `http://sun.stanford.edu/~rmunk/PROPACK/` while he was at
> Stanford University. This documentation preserves the material from that
> original homepage and keeps it alongside the modern library in this
> repository. See [history and provenance](history.md).

## The modern library

The actively maintained Fortran 77 / C / C++ implementation lives in this
repository. See the top-level [README](../README.md) for build instructions
(CMake), the header-only C++17 API, and benchmarks.

Two SVD algorithm variants are provided:

- **`Xlansvd`** — basic Lanczos bidiagonalization (memory grows with the number
  of iterations, but typically uses fewer matrix-vector products).
- **`Xlansvd_irl`** — implicitly restarted variant, which computes a given
  number of singular triplets in a *fixed* amount of memory.

where `X` is `d`, `s`, `z`, or `c` for double, single, complex-double, and
complex-single precision respectively.

## Documentation

- [Technical report and talks](papers/index.md) — the 1998 report that
  introduced BPRO, plus two conference presentations, with abstracts.
- [Algorithm and routine reference](routines/index.md) — the four core
  routines: `lansvd`, `lanbpro`, `laneig`, `lanpro`.
- [Fortran package guide (v2.1)](guides/fortran-v2.1.md) — overview and
  installation notes for the classic Fortran 77 distribution.
- [MATLAB toolbox guide](guides/matlab.md) — overview of the original MATLAB
  implementation and its computational routines.
- [References](references.md) — bibliography.
- [History and provenance](history.md) — where this material came from and what
  was corrected.
- [Legacy archive](legacy/index.md) — verbatim copies of the files from the
  original Stanford homepage.

## The algorithm in brief

For an `m`-by-`n` matrix `A`, `k` steps of Lanczos bidiagonalization produce

```
A V_k = U_k B_k + r_k e_k^T
```

with `B_k` lower bidiagonal and `U_k`, `V_k` having (semi-)orthonormal columns.
The singular values of `B_k` approximate those of `A`. Maintaining full
orthogonality of the Lanczos vectors is expensive; **partial
reorthogonalization** monitors the loss of orthogonality through simple scalar
recurrences — without forming inner products — and reorthogonalizes only when
semi-orthogonality is about to be lost. This is mathematically equivalent to
Simon's symmetric Lanczos PRO, but works directly on the bidiagonalization of
`A`, using about half the storage and a factor of 3–4 less work than applying
PRO to an equivalent symmetric system.

## Citation

If you use PROPACK, please cite the technical report:

> R. M. Larsen. *Lanczos bidiagonalization with partial reorthogonalization.*
> Department of Computer Science, Århus University, Technical Report DAIMI
> PB-537, October 1998.
> DOI: [10.7146/dpb.v27i537.7070](https://doi.org/10.7146/dpb.v27i537.7070)

A machine-readable [`CITATION.cff`](../CITATION.cff) is provided at the
repository root.

## License

- **Software:** BSD 3-Clause — see [`LICENSE`](../LICENSE).
- **Documentation** in this `docs/` tree: Creative Commons Attribution 4.0
  (CC BY 4.0) — see [`LICENSE-docs.md`](LICENSE-docs.md) — except for the
  verbatim files under [`legacy/`](legacy/index.md), which retain their original
  notices.

## Contact

Rasmus Munk Larsen — rmlarsen@gmail.com
