# Algorithm and routine reference

PROPACK is built around four core routines. The two Lanczos drivers (`lanbpro`,
`lanpro`) implement the iterative bidiagonalization / tridiagonalization with
partial reorthogonalization; the two solvers (`lansvd`, `laneig`) build on them
to compute singular triplets and eigenpairs.

| Routine | Purpose |
|---------|---------|
| [`lansvd`](lansvd.md)   | A few singular values and singular vectors of a general matrix (SVD), via Lanczos **bi**diagonalization with partial reorthogonalization (BPRO). |
| [`lanbpro`](lanbpro.md) | Lanczos **bi**diagonalization with partial reorthogonalization — the engine underlying `lansvd`. |
| [`laneig`](laneig.md)   | A few eigenvalues and eigenvectors of a real symmetric matrix, via the symmetric Lanczos algorithm with partial reorthogonalization (PRO). |
| [`lanpro`](lanpro.md)   | Lanczos **tri**diagonalization with partial reorthogonalization — the engine underlying `laneig`. |

The pages below reproduce the on-line help for the original MATLAB routines. In
the modern Fortran library the SVD drivers are named `Xlansvd` and `Xlansvd_irl`
(`X` ∈ {`d`, `s`, `z`, `c`}); see the [Fortran guide](../guides/fortran-v2.1.md)
and the top-level [README](../../README.md).

## Common reorthogonalization options

All four routines accept an `OPTIONS` structure controlling reorthogonalization.
Exact defaults vary slightly by routine; see each page.

| Field | Meaning |
|-------|---------|
| `delta` | Desired level of orthogonality among the Lanczos vectors (default `sqrt(eps/K)`). |
| `eta`   | Level of orthogonality enforced after reorthogonalization. |
| `cgs`   | Reorthogonalization method: `0` = iterated modified Gram-Schmidt (default), `1` = iterated classical Gram-Schmidt. |
| `elr`   | If `1` (default), extended local reorthogonalization is enforced. |

See [references](../references.md) for the underlying theory.
