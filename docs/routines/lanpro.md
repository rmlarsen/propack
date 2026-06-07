---
title: "lanpro — Lanczos tridiagonalization (PROPACK)"
description: "Lanczos tridiagonalization with partial reorthogonalization, the engine underlying laneig."
---

# `lanpro` — Lanczos tridiagonalization with partial reorthogonalization

`LANPRO` computes the Lanczos tridiagonalization of a real symmetric matrix using
the symmetric Lanczos algorithm with partial reorthogonalization (PRO). It is the
engine underlying [`laneig`](laneig.md).

## Syntax

```matlab
[Q_K,T_K,R,ANORM,IERR,WORK] = lanpro(A,K,R0,OPTIONS,Q_old,T_old)
[Q_K,T_K,R,ANORM,IERR,WORK] = lanpro('Afun',N,K,R0,OPTIONS,Q_old,T_old)
```

Computes `K` steps of the Lanczos algorithm with starting vector `R0`, returning
the `K`-by-`K` tridiagonal `T_K`, the `N`-by-`K` matrix `Q_K` with
semiorthonormal columns, and the residual vector `R` such that

```
A*Q_K = Q_K*T_K + R
```

Partial reorthogonalization keeps the columns of `Q_K` semiorthogonal:

```
max(max(abs(eye(K) - Q_K'*Q_K))) <= OPTIONS.delta
```

The first input is either a real symmetric matrix, a struct with components `A.L`
and `A.U`, or a string naming an M-file that applies a linear operator to the
columns of a given matrix. In the latter case, the second input must be `N`, the
order of the problem. If `A` is a struct with components `A.L` and `A.U` such
that `L*U = (A - sigma*I)`, a shift-and-invert Lanczos iteration is performed.

## Options

| Field | Meaning | Default |
|-------|---------|---------|
| `OPTIONS.delta` | Desired level of orthogonality | `sqrt(eps/K)` |
| `OPTIONS.eta`   | Level of orthogonality after reorthogonalization | `eps^(3/4)/sqrt(K)` |
| `OPTIONS.cgs`   | `0` = iterated modified Gram-Schmidt, `1` = iterated classical Gram-Schmidt | `0` |
| `OPTIONS.elr`   | If `1`, extended local reorthogonalization is enforced | `1` |
| `OPTIONS.Y`     | The Lanczos vectors are reorthogonalized against the columns of `OPTIONS.Y` | — |

## Extending a factorization

If `R0`, `Q_old`, and `T_old` are provided, they must contain a partial Lanczos
tridiagonalization of `A` of the form `A*Q_old = Q_old*T_old + R0`. The
factorization is extended to dimension `K`-by-`K` by continuing the iteration
with `R0` as the starting vector.

## Outputs

On exit, `ANORM` contains an approximation to `||A||_2`, and:

| `IERR` | Meaning |
|--------|---------|
| `= 0`  | `K` steps were performed successfully. |
| `> 0`  | `K` steps were performed, but the algorithm switched to full reorthogonalization after `IERR` steps. |
| `< 0`  | Iteration terminated after `-IERR` steps because an invariant subspace was found and 3 deflation attempts were unsuccessful. |

`WORK(1)` contains the number of reorthogonalizations performed, and `WORK(2)`
the number of inner products performed in the reorthogonalizations.

## See also

[`laneig`](laneig.md), and the computational routines `reorth`, `compute_int`.

## References

See [references](../references.md). Key sources: R. M. Larsen (Århus, 1998);
G. H. Golub & C. F. Van Loan, *Matrix Computations*, 3rd ed., Johns Hopkins,
1996 (Ch. 9); B. N. Parlett (1980); H. D. Simon, *Math. Comp.* 42 (1984),
115–142.
