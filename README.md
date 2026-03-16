# PROPACK

A Fortran library for computing the partial singular value decomposition (SVD)
of large, sparse, and/or structured matrices via Lanczos bidiagonalization with
partial reorthogonalization.

Two algorithm variants are provided:

- **`Xlansvd`** -- basic Lanczos bidiagonalization (memory grows with
  iteration count, but typically uses fewer matrix-vector products).
- **`Xlansvd_irl`** -- implicitly restarted variant (fixed memory usage).

where `X` is `d`, `s`, `z`, or `c` for double, single, complex double, and
complex single precision, respectively.

## Quick start

```sh
cmake -B build
cmake --build build
ctest --test-dir build
```

## CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `PROPACK_BUILD_EXAMPLES` | `ON` | Build example programs |
| `PROPACK_BUILD_BENCHMARKS` | `OFF` | Build benchmark programs |
| `PROPACK_BUILD_TESTING` | `ON` | Build CTest targets for examples |
| `PROPACK_USE_OPENMP` | `OFF` | Enable OpenMP parallelization |
| `PROPACK_BUILD_SINGLE` | `ON` | Build single-precision (s) library |
| `PROPACK_BUILD_DOUBLE` | `ON` | Build double-precision (d) library |
| `PROPACK_BUILD_COMPLEX8` | `ON` | Build complex single-precision (c) library |
| `PROPACK_BUILD_COMPLEX16` | `ON` | Build complex double-precision (z) library |
| `PROPACK_BUILD_CXX_API` | `ON` | Build C++17 API wrappers (`propack.hpp`) |

Example:

```sh
cmake -B build -DPROPACK_USE_OPENMP=ON -DPROPACK_BUILD_BENCHMARKS=ON
```

## BLAS/LAPACK

CMake auto-detects BLAS and LAPACK. OpenBLAS is preferred when available;
otherwise any BLAS/LAPACK found by CMake's `find_package(LAPACK)` is used.

Override with:

```sh
cmake -B build -DBLA_VENDOR=Intel10_64lp_seq
```

See the [CMake FindLAPACK documentation](https://cmake.org/cmake/help/latest/module/FindLAPACK.html)
for the full list of `BLA_VENDOR` values.

## Supported compilers

The build system provides tuned default flags for:

- GNU (`gfortran`)
- Intel classic (`ifort`)
- Intel LLVM (`ifx`)
- LLVM Flang (`flang-new` / `flang`)
- NVIDIA HPC SDK (`nvfortran`)
- NAG (`nagfor`)

Other Fortran compilers should work but may need manual flag tuning via
`-DCMAKE_Fortran_FLAGS=...`.

## Precision variants

Four static libraries are built (one per precision):

| Library | Precision | Alias target |
|---------|-----------|--------------|
| `libdpropack.a` | double (`real*8`) | `PROPACK::dpropack` |
| `libspropack.a` | single (`real*4`) | `PROPACK::spropack` |
| `libzpropack.a` | complex double (`complex*16`) | `PROPACK::zpropack` |
| `libcpropack.a` | complex single (`complex*8`) | `PROPACK::cpropack` |

A shared `propack_common` library (integer utilities, statistics printing) is
linked automatically.

## Example programs

Each precision directory contains an `Examples/` subdirectory with:

- **`example`** -- compute a partial SVD using `Xlansvd`.
- **`example_irl`** -- compute a partial SVD using `Xlansvd_irl`.

The example programs read a sparse matrix in Harwell-Boeing format from a file.
Test matrices are provided: `illc1850.rra` (real precisions) and `mhd1280b.cua`
(complex precisions). The programs also support diagonal, coordinate, and dense
matrix formats.

## Testing

```sh
ctest --test-dir build
```

Tests use CTest fixture dependencies to ensure correct ordering: example
programs run first to produce output, then the `compare` program verifies
results against reference output in `<precision>/Examples/Output/`.

## Benchmarks

Enable benchmarks with:

```sh
cmake -B build -DPROPACK_BUILD_BENCHMARKS=ON
cmake --build build
```

Run individually via CTest:

```sh
ctest --test-dir build -L benchmark
```

Or use the helper script for formatted Markdown table output:

```sh
./run_benchmarks.sh build
```

### SuiteSparse matrix benchmarks

Download the 20 SuiteSparse test matrices used in the propack++ paper:

```sh
./paper/download_matrices.sh paper/matrices
```

Then run the Matrix Market benchmark (requires Eigen3):

```sh
cmake -B build -DPROPACK_BUILD_BENCHMARKS=ON
cmake --build build
./build/propack_bench_mtx paper/matrices 10 50
```

This benchmarks both `lansvd` and `lansvd_irl` with CGS and MGS on all
downloaded matrices, outputting JSON results to stdout.

## Using PROPACK in your project

After installing PROPACK (`cmake --install build`), use it from another CMake
project:

```cmake
find_package(PROPACK REQUIRED)
target_link_libraries(mytarget PRIVATE PROPACK::dpropack)
```

## API

The user supplies an `APROD` subroutine that computes matrix-vector products.
Its interface is:

```fortran
SUBROUTINE APROD(TRANSA, M, N, X, Y, DPARM, IPARM)
  CHARACTER*1 TRANSA
  INTEGER     M, N, IPARM(*)
  <type>      X(*), Y(*), DPARM(*)
```

- `TRANSA = 'N'`: compute `Y = A * X`
- `TRANSA = 'T'` (real) or `'C'` (complex): compute `Y = A^T * X` or `Y = A^H * X`
- `DPARM` and `IPARM` pass user data to `APROD` without common blocks.

See the source files `<precision>/<X>lansvd.F` and `<precision>/<X>lansvd_irl.F`
for full parameter documentation.

### Option arrays

Algorithm behavior is controlled via `DOPTION` and `IOPTION` arrays:

| Entry | Default | Description |
|-------|---------|-------------|
| `doption(1)` | `sqrt(eps)` | `delta` — level of orthogonality to maintain |
| `doption(2)` | `eps^(3/4)` | `eta` — reorthogonalization purge threshold |
| `doption(3)` | `0` (auto) | `anorm` — estimate of ‖A‖ |
| `doption(4)` | `1e-3` | min relative gap between shifts (IRL only) |
| `doption(5)` | `0` (auto) | `eps1_aprod` — rounding error for A*v (auto: `5*sqrt(n)*eps`) |
| `doption(6)` | `0` (auto) | `eps1_atprod` — rounding error for A'*u (auto: `5*sqrt(m)*eps`) |
| `ioption(1)` | `0` | `MGS` — use modified Gram-Schmidt (0=CGS, 1=MGS) |
| `ioption(2)` | `1` | `ELR` — extended local reorthogonalization |

**Important:** `doption` must have at least 6 elements, all initialized.
Set entries 5-6 to 0 for automatic (recommended) values.

### Error codes

| `INFO` value | Meaning |
|-------------|---------|
| `0` | Success — all requested singular triplets converged |
| `> 0` | An invariant subspace of dimension `INFO` was found |
| `< 0, > -1000` | Only `|INFO|` singular triplets converged |
| `-1000` | LAPACK bidiagonal SVD failed |
| `-2000` | Ghost singular value detected — partial reorthogonalization failed |

## C++ API

When `PROPACK_BUILD_CXX_API=ON` (the default), a header-only C++17 interface is
installed as `propack.hpp`. It provides `propack::lansvd<T>()` and
`propack::lansvd_irl<T>()` for `T` = `double`, `float`, `std::complex<double>`,
or `std::complex<float>`. Workspace is managed automatically. Requires a C++17
compiler.

```cpp
#include <propack.hpp>
#include <iostream>

int main() {
    const int m = 100, n = 80, k = 5, kmax = 20;
    // Diagonal matrix with entries 1, 2, ..., min(m,n)
    auto result = propack::lansvd<double>(m, n, k, kmax,
        [](char trans, int m, int n, const double* x, double* y) {
            int len = (trans == 'N') ? m : n;
            int mn  = (m < n) ? m : n;
            for (int i = 0; i < len; ++i) y[i] = 0.0;
            for (int i = 0; i < mn; ++i)
                y[i] = (i + 1.0) * x[i];  // diag(1,2,...,min(m,n))
        });

    std::cout << "info = " << result.info << "\n";
    for (int i = 0; i < k; ++i)
        std::cout << "sigma[" << i << "] = " << result.sigma[i]
                  << " +/- " << result.bnd[i] << "\n";
}
```

Link against the appropriate PROPACK library and a Fortran runtime:

```cmake
find_package(PROPACK REQUIRED)
target_link_libraries(myapp PRIVATE PROPACK::dpropack)
```

## Contact

Questions and comments about PROPACK are welcome and should be directed to:

Rasmus Munk Larsen
E-mail: rmlarsen@gmail.com

(C) Rasmus Munk Larsen, Mountain View, October 2025.
