---
title: "lansvd — partial SVD (PROPACK)"
description: "Compute a few singular values and vectors via Lanczos bidiagonalization with partial reorthogonalization."
---

# `lansvd` — singular values and vectors via Lanczos bidiagonalization

`LANSVD` computes singular triplets `(u, v, sigma)` such that `A*u = sigma*v` and
`A'*v = sigma*u`. Only a few singular values and singular vectors are computed,
using the Lanczos bidiagonalization algorithm with partial reorthogonalization
(BPRO).

## Syntax

```matlab
S = lansvd(A)
S = lansvd('Afun','Atransfun',M,N)

[U,S,V] = lansvd(A,K,'L',...)   % K largest singular values
[U,S,V] = lansvd(A,K,'S',...)   % K smallest singular values

[U,S,V] = lansvd(A,K,SIGMA,OPTIONS)
[U,S,V] = lansvd('Afun','Atransfun',M,N,K,SIGMA,OPTIONS)
```

The first input argument is either a matrix or a string naming an M-file that
applies a linear operator to the columns of a given matrix. In the latter case,
the second input must name an M-file that applies the *transpose* of the same
operator, and the third and fourth arguments must be `M` and `N`, the dimensions
of the problem. `K` is the number of singular values desired and `SIGMA` is `'L'`
(largest) or `'S'` (smallest).

## Options

| Field | Parameter | Default |
|-------|-----------|---------|
| `OPTIONS.tol`    | Convergence tolerance | `16*eps` |
| `OPTIONS.lanmax` | Dimension of the Lanczos basis | — |
| `OPTIONS.p0`     | Starting vector for the Lanczos iteration | `rand(n,1)-0.5` |
| `OPTIONS.delta`  | Level of orthogonality among the Lanczos vectors | `sqrt(eps/K)` |
| `OPTIONS.eta`    | Level of orthogonality after reorthogonalization | `10*eps^(3/4)` |
| `OPTIONS.cgs`    | `0` = iterated modified Gram-Schmidt, `1` = iterated classical Gram-Schmidt | `0` |
| `OPTIONS.elr`    | If `1`, extended local reorthogonalization is enforced | `1` |

## See also

[`lanbpro`](lanbpro.md), and the MATLAB built-ins `svds`, `svd`.

## References

See [references](../references.md). Key sources: R. M. Larsen, *Lanczos
bidiagonalization with partial reorthogonalization* (Århus University, 1998);
H. D. Simon, "The Lanczos algorithm with partial reorthogonalization",
*Math. Comp.* 42 (1984), 115–142; B. N. Parlett, *The Symmetric Eigenvalue
Problem*, Prentice-Hall, 1980.
