---
title: "lanbpro — Lanczos bidiagonalization (PROPACK)"
description: "Lanczos bidiagonalization with partial reorthogonalization, the engine underlying lansvd."
---

# `lanbpro` — Lanczos bidiagonalization with partial reorthogonalization

`LANBPRO` computes the Lanczos bidiagonalization of a rectangular matrix using
partial reorthogonalization. It is the engine underlying [`lansvd`](lansvd.md).

## Syntax

```matlab
[U_k,B_k,V_k,R,ierr,work] = lanbpro(A,K,R0,OPTIONS,U_old,B_old,V_old)
[U_k,B_k,V_k,R,ierr,work] = lanbpro('Afun','Atransfun',M,N,K,R0, ...
                                    OPTIONS,U_old,B_old,V_old)

B_k = lanbpro(...)   % returns the bidiagonal matrix only
```

Computes `K` steps of the Lanczos bidiagonalization algorithm with partial
reorthogonalization (BPRO) on the `M`-by-`N` matrix `A` with `M`-by-1 starting
vector `R0`, producing a lower bidiagonal `K`-by-`K` matrix `B_k`, an `N`-by-`K`
matrix `V_k`, an `M`-by-`K` matrix `U_k`, and an `M`-by-1 vector `R` such that

```
A*V_k = U_k*B_k + R
```

Partial reorthogonalization keeps the columns of `V_k` and `U_k` semiorthogonal:

```
max(max(abs(eye(K) - V_k'*V_k))) <= OPTIONS.delta
max(max(abs(eye(K) - U_k'*U_k))) <= OPTIONS.delta
```

The first input is either a real matrix or a string naming an M-file that applies
a linear operator to the columns of a given matrix. In the latter case, the
second input must name an M-file applying the *transpose* of the same operator,
and the third and fourth arguments must be `M` and `N`.

## Options

| Field | Meaning | Default |
|-------|---------|---------|
| `OPTIONS.delta`    | Desired level of orthogonality | `sqrt(eps/K)` |
| `OPTIONS.eta`      | Level of orthogonality after reorthogonalization | `eps^(3/4)/sqrt(K)` |
| `OPTIONS.cgs`      | `0` = iterated modified Gram-Schmidt, `1` = iterated classical Gram-Schmidt | `0` |
| `OPTIONS.elr`      | If `1`, extended local reorthogonalization is enforced | `1` |
| `OPTIONS.onesided` | `0` = keep both `U` and `V` semiorthogonal; `1` = reorthogonalize only `U`; `-1` = reorthogonalize only `V` | `0` |
| `OPTIONS.waitbar`  | If set, the progress of the algorithm is displayed graphically | — |

## Extending a factorization

If `R0`, `U_old`, `B_old`, and `V_old` are all provided, they must contain a
partial Lanczos bidiagonalization of `A` of the form `A*V_old = U_old*B_old + R0`.
The factorization is then extended to dimension `K`-by-`K` by continuing the
Lanczos bidiagonalization with `R0` as the starting vector.

## Output `work`

The output array `work` contains information about the reorthogonalization work:

```
work = [ RU  PU ]
       [ RV  PV ]
```

where `RU`/`RV` are the number of reorthogonalizations of `U`/`V`, and `PU`/`PV`
are the number of inner products used in reorthogonalizing `U`/`V`.

## References

See [references](../references.md). Key sources: R. M. Larsen (Århus, 1998);
G. H. Golub & C. F. Van Loan, *Matrix Computations*, 3rd ed., Johns Hopkins,
1996 (§9.3.4); B. N. Parlett (1980); H. D. Simon, *Math. Comp.* 42 (1984),
115–142.
