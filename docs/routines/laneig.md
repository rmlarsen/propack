---
title: "laneig — symmetric eigenproblem (PROPACK)"
description: "A few eigenvalues and eigenvectors of a real symmetric matrix via Lanczos with partial reorthogonalization."
---

# `laneig` — eigenvalues and eigenvectors via the symmetric Lanczos algorithm

`LANEIG` solves the eigenvalue problem `A*v = lambda*v` when `A` is real and
symmetric, using the Lanczos algorithm with partial reorthogonalization (PRO).
It builds on [`lanpro`](lanpro.md).

## Syntax

```matlab
[V,D] = laneig(A)
[V,D] = laneig('Afun',N)

[V,D,ERR] = laneig(A,K,SIGMA,OPTIONS)
[V,D,ERR] = laneig('Afun',N,K,SIGMA,OPTIONS)
```

The first input is either a real symmetric matrix or a string naming an M-file
that applies a linear operator to the columns of a given matrix. In the latter
case, the second input must be `N`, the order of the problem. On exit, `ERR`
contains the computed error bounds. `K` is the number of eigenvalues desired.

`SIGMA` is a numerical shift or a two-letter string specifying which part of the
spectrum to compute:

| `SIGMA` | Eigenvalues computed |
|---------|----------------------|
| `'AL'`  | Algebraically largest |
| `'AS'`  | Algebraically smallest |
| `'LM'`  | Largest magnitude (default) |
| `'SM'`  | Smallest magnitude (does not work when `A` is an M-file) |
| `'BE'`  | Both ends — `k/2` eigenvalues from each end of the spectrum (one extra from the high end if `k` is odd) |

## Options

| Field | Parameter | Default |
|-------|-----------|---------|
| `OPTIONS.tol`    | Convergence tolerance | `16*eps` |
| `OPTIONS.lanmax` | Dimension of the Lanczos basis | — |
| `OPTIONS.v0`     | Starting vector for the Lanczos iteration | `rand(n,1)-0.5` |
| `OPTIONS.delta`  | Level of orthogonality among the Lanczos vectors | `sqrt(eps/K)` |
| `OPTIONS.eta`    | Level of orthogonality after reorthogonalization | `10*eps^(3/4)` |
| `OPTIONS.cgs`    | `0` = iterated modified Gram-Schmidt, `1` = iterated classical Gram-Schmidt | `0` |
| `OPTIONS.elr`    | If `1`, extended local reorthogonalization is enforced | `1` |

## See also

[`lanpro`](lanpro.md), and the MATLAB built-ins `eigs`, `eig`.

## References

See [references](../references.md). Key sources: R. M. Larsen (Århus, 1998);
B. N. Parlett, *The Symmetric Eigenvalue Problem*, Prentice-Hall, 1980;
H. D. Simon, *Math. Comp.* 42 (1984), 115–142.
