// bench_propack_cxx.cpp — Google Benchmark suite for PROPACK C++ API
//
// Mirrors the Fortran benchmark_gen.F structure using two matrix types
// (random sparse CCS and diagonal with power-law decay) at various sizes.
// Double precision only.

#include <benchmark/benchmark.h>
#include <propack.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// ============================================================================
// Park-Miller LCG (matching Fortran benchmark_gen.F lcg_next)
// ============================================================================

static int32_t lcg_next(int32_t seed) {
    int64_t s = seed;
    s = (s * 16807) % INT64_C(2147483647);
    return static_cast<int32_t>(s);
}

// ============================================================================
// Sparse matrix in CCS (Compressed Column Storage) format
// ============================================================================

struct SparseMatrix {
    std::vector<int> colptr;
    std::vector<int> rowind;
    std::vector<double> values;
    int m, n;
};

// Generate random sparse CCS matrix matching Fortran gen_random_sparse.
// Row indices are 0-based (Fortran uses 1-based).
static SparseMatrix gen_random_sparse(int m, int n, int nnz_per_col, int seed) {
    SparseMatrix mat;
    mat.m = m;
    mat.n = n;
    int nnz = n * nnz_per_col;
    mat.colptr.resize(n + 1);
    mat.rowind.resize(nnz);
    mat.values.resize(nnz);

    int32_t iseed = seed;
    int idx = 0;
    for (int j = 0; j < n; ++j) {
        mat.colptr[j] = idx;
        for (int kk = 0; kk < nnz_per_col; ++kk) {
            iseed = lcg_next(iseed);
            mat.rowind[idx] = std::abs(iseed) % m;
            iseed = lcg_next(iseed);
            mat.values[idx] =
                2.0 * static_cast<double>(iseed) / 2147483647.0 - 1.0;
            ++idx;
        }
    }
    mat.colptr[n] = idx;
    return mat;
}

// ============================================================================
// Diagonal vector with power-law singular value decay
// diag(i) = (min(m,n) / (i+1))^alpha   (0-indexed)
// ============================================================================

static std::vector<double> gen_diagonal(int m, int n, double alpha) {
    int mn = std::min(m, n);
    std::vector<double> diag(mn);
    double dmn = static_cast<double>(mn);
    for (int i = 0; i < mn; ++i) {
        diag[i] = std::pow(dmn / static_cast<double>(i + 1), alpha);
    }
    return diag;
}

// ============================================================================
// APROD implementations — return lambda capturing data via shared_ptr
// ============================================================================

template <typename Scalar>
propack::AprodFn<Scalar>
make_sparse_aprod(std::shared_ptr<SparseMatrix> mat) {
    return [mat](char trans, int m, int n, const Scalar* x, Scalar* y) {
        if (trans == 'N' || trans == 'n') {
            // y = A * x
            std::fill(y, y + m, Scalar(0));
            for (int j = 0; j < n; ++j) {
                for (int kk = mat->colptr[j]; kk < mat->colptr[j + 1]; ++kk) {
                    y[mat->rowind[kk]] +=
                        static_cast<Scalar>(mat->values[kk]) * x[j];
                }
            }
        } else {
            // y = A^T * x
            std::fill(y, y + n, Scalar(0));
            for (int j = 0; j < n; ++j) {
                for (int kk = mat->colptr[j]; kk < mat->colptr[j + 1]; ++kk) {
                    y[j] += static_cast<Scalar>(mat->values[kk]) *
                            x[mat->rowind[kk]];
                }
            }
        }
    };
}

template <typename Scalar>
propack::AprodFn<Scalar>
make_diagonal_aprod(std::shared_ptr<std::vector<double>> diag) {
    return [diag](char trans, int m, int n, const Scalar* x, Scalar* y) {
        int mn = static_cast<int>(diag->size());
        int outlen = (trans == 'N' || trans == 'n') ? m : n;
        std::fill(y, y + outlen, Scalar(0));
        for (int i = 0; i < mn; ++i) {
            y[i] = static_cast<Scalar>((*diag)[i]) * x[i];
        }
    };
}

// ============================================================================
// Benchmark registration
// ============================================================================

static void register_all() {
    // Common options matching Fortran: extended_local_reorth=true, CGS reorth
    propack::Options opts;
    opts.extended_local_reorth = true;

    propack::IRLOptions irl_opts;
    irl_opts.extended_local_reorth = true;

    constexpr double tol = 1e-12;
    constexpr int kmax_limit = 300;  // matches Fortran gkmax
    constexpr int maxiter = 1000;

    // Compute IRL parameters matching Fortran benchmark_gen.F
    auto irl_params = [](int k, int m, int n) -> std::pair<int, int> {
        int dim = std::min(2 * k + 10, std::min(m, n));
        int p = std::min(dim - k, k);
        return {dim, p};
    };

    // -----------------------------------------------------------------
    // Sparse configurations (nnz_per_col=20, seed=12345)
    // -----------------------------------------------------------------
    struct Config {
        int m, n, k;
        const char* label;
    };
    Config sparse_cfgs[] = {
        {1000, 1000,  50, "1000x1000_k50"},
        {2000, 2000,  50, "2000x2000_k50"},
        {3000, 3000,  50, "3000x3000_k50"},
        {3000, 3000, 100, "3000x3000_k100"},
        {5000, 5000,  50, "5000x5000_k50"},
        {5000, 5000, 100, "5000x5000_k100"},
        {5000, 2000,  50, "5000x2000_k50"},
    };

    for (const auto& c : sparse_cfgs) {
        auto mat = std::make_shared<SparseMatrix>(
            gen_random_sparse(c.m, c.n, 20, 12345));
        int m = c.m, n = c.n, k = c.k;
        int kmax = std::min(kmax_limit, std::min(m, n));
        auto [dim, p] = irl_params(k, m, n);

        benchmark::RegisterBenchmark(
            (std::string("BM_lansvd/sparse_") + c.label).c_str(),
            [=](benchmark::State& state) {
                auto aprod = make_sparse_aprod<double>(mat);
                for (auto _ : state) {
                    auto r = propack::lansvd<double>(
                        m, n, k, kmax, aprod, false, false, tol, opts);
                    benchmark::DoNotOptimize(r);
                }
            });

        benchmark::RegisterBenchmark(
            (std::string("BM_lansvd_irl/sparse_") + c.label).c_str(),
            [=](benchmark::State& state) {
                auto aprod = make_sparse_aprod<double>(mat);
                for (auto _ : state) {
                    auto r = propack::lansvd_irl<double>(
                        'L', m, n, k, dim, p, maxiter, aprod,
                        false, false, tol, irl_opts);
                    benchmark::DoNotOptimize(r);
                }
            });
    }

    // -----------------------------------------------------------------
    // Diagonal configurations (alpha=1.0)
    // -----------------------------------------------------------------
    Config diag_cfgs[] = {
        {2000, 2000,  50, "2000x2000_k50"},
        {2000, 2000, 100, "2000x2000_k100"},
        {5000, 5000,  50, "5000x5000_k50"},
        {5000, 5000, 100, "5000x5000_k100"},
    };

    for (const auto& c : diag_cfgs) {
        auto diag = std::make_shared<std::vector<double>>(
            gen_diagonal(c.m, c.n, 1.0));
        int m = c.m, n = c.n, k = c.k;
        int kmax = std::min(kmax_limit, std::min(m, n));
        auto [dim, p] = irl_params(k, m, n);

        benchmark::RegisterBenchmark(
            (std::string("BM_lansvd/diag_") + c.label).c_str(),
            [=](benchmark::State& state) {
                auto aprod = make_diagonal_aprod<double>(diag);
                for (auto _ : state) {
                    auto r = propack::lansvd<double>(
                        m, n, k, kmax, aprod, false, false, tol, opts);
                    benchmark::DoNotOptimize(r);
                }
            });

        benchmark::RegisterBenchmark(
            (std::string("BM_lansvd_irl/diag_") + c.label).c_str(),
            [=](benchmark::State& state) {
                auto aprod = make_diagonal_aprod<double>(diag);
                for (auto _ : state) {
                    auto r = propack::lansvd_irl<double>(
                        'L', m, n, k, dim, p, maxiter, aprod,
                        false, false, tol, irl_opts);
                    benchmark::DoNotOptimize(r);
                }
            });
    }
}

// Register benchmarks at static-init time (before benchmark_main's main()).
namespace {
struct Registrar {
    Registrar() { register_all(); }
} registrar;
}
