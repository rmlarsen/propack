// Unit tests for the PROPACK C++ API (propack.hpp).
// Uses Google Test with typed tests to cover all four scalar types.

#include <propack.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

// ============================================================================
// Helpers
// ============================================================================

template<typename T> struct RealOf              { using type = T; };
template<typename T> struct RealOf<std::complex<T>> { using type = T; };

// Tolerances per precision.
template<typename T> struct Tol;
template<> struct Tol<double>                { static constexpr double v = 1e-6;  };
template<> struct Tol<float>                 { static constexpr float  v = 5e-4f; };
template<> struct Tol<std::complex<double>>  { static constexpr double v = 1e-6;  };
template<> struct Tol<std::complex<float>>   { static constexpr float  v = 5e-4f; };

// Return an APROD lambda for diag(1,2,...,min(m,n)).
template<typename Scalar>
propack::AprodFn<Scalar> diagonal_aprod() {
    return [](char trans, int m, int n, const Scalar* x, Scalar* y) {
        int len = (trans == 'N' || trans == 'n') ? m : n;
        int mn  = std::min(m, n);
        for (int i = 0; i < len; ++i) y[i] = Scalar(0);
        for (int i = 0; i < mn; ++i)
            y[i] = static_cast<typename RealOf<Scalar>::type>(i + 1) * x[i];
    };
}

// Compute max |col_i^H col_j - delta_{ij}| for a column-major matrix.
template<typename Scalar>
double orthogonality_error(const Scalar* data, int rows, int cols) {
    using Real = typename RealOf<Scalar>::type;
    double max_err = 0.0;
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i <= j; ++i) {
            // dot(col_i, col_j)
            std::complex<double> dot(0.0);
            for (int r = 0; r < rows; ++r) {
                auto a = static_cast<std::complex<double>>(data[i * rows + r]);
                auto b = static_cast<std::complex<double>>(data[j * rows + r]);
                dot += std::conj(a) * b;
            }
            double expected = (i == j) ? 1.0 : 0.0;
            double err = std::abs(dot - std::complex<double>(expected));
            max_err = std::max(max_err, err);
        }
    }
    return max_err;
}

// ============================================================================
// Typed test fixture — parametrized on Scalar
// ============================================================================

template<typename Scalar>
class LansvdTest : public ::testing::Test {};

using ScalarTypes = ::testing::Types<double, float,
                                     std::complex<double>,
                                     std::complex<float>>;
TYPED_TEST_SUITE(LansvdTest, ScalarTypes);

// ============================================================================
// lansvd: diagonal matrix, verify singular values
// ============================================================================

TYPED_TEST(LansvdTest, DiagonalSingularValues) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, k = 5, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>());

    EXPECT_EQ(result.info, 0);
    EXPECT_EQ(static_cast<int>(result.sigma.size()), k);
    EXPECT_EQ(static_cast<int>(result.bnd.size()), k);

    for (int i = 0; i < k; ++i) {
        Real expected = static_cast<Real>(n - i);  // 80, 79, 78, 77, 76
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v))
            << "sigma[" << i << "]";
    }
}

// ============================================================================
// lansvd: verify error bounds are small
// ============================================================================

TYPED_TEST(LansvdTest, ErrorBoundsSmall) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, k = 5, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>());

    ASSERT_EQ(result.info, 0);
    for (int i = 0; i < k; ++i) {
        EXPECT_LT(static_cast<double>(result.bnd[i]), 1e-3)
            << "bnd[" << i << "] too large";
    }
}

// ============================================================================
// lansvd: U has correct dimensions, columns are orthonormal
// ============================================================================

TYPED_TEST(LansvdTest, LeftVectorsOrthonormal) {
    using Scalar = TypeParam;

    const int m = 100, n = 80, k = 5, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>(),
                                          /*compute_u=*/true, /*compute_v=*/false);

    ASSERT_EQ(result.info, 0);
    EXPECT_EQ(static_cast<int>(result.U.size()), m * k);
    EXPECT_TRUE(result.V.empty());

    double orth_err = orthogonality_error(result.U.data(), m, k);
    EXPECT_LT(orth_err, static_cast<double>(Tol<Scalar>::v))
        << "U columns not orthonormal";
}

// ============================================================================
// lansvd: V has correct dimensions, columns are orthonormal
// ============================================================================

TYPED_TEST(LansvdTest, RightVectorsOrthonormal) {
    using Scalar = TypeParam;

    const int m = 100, n = 80, k = 5, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>(),
                                          /*compute_u=*/false, /*compute_v=*/true);

    ASSERT_EQ(result.info, 0);
    EXPECT_TRUE(result.U.empty());
    EXPECT_EQ(static_cast<int>(result.V.size()), n * k);

    double orth_err = orthogonality_error(result.V.data(), n, k);
    EXPECT_LT(orth_err, static_cast<double>(Tol<Scalar>::v))
        << "V columns not orthonormal";
}

// ============================================================================
// lansvd: no vectors mode — U and V empty
// ============================================================================

TYPED_TEST(LansvdTest, NoVectors) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, k = 5, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>(),
                                          /*compute_u=*/false, /*compute_v=*/false);

    EXPECT_EQ(result.info, 0);
    EXPECT_TRUE(result.U.empty());
    EXPECT_TRUE(result.V.empty());
    // Singular values should still be correct.
    for (int i = 0; i < k; ++i) {
        Real expected = static_cast<Real>(n - i);
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v));
    }
}

// ============================================================================
// lansvd: tall matrix (m > n)
// ============================================================================

TYPED_TEST(LansvdTest, TallMatrix) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 120, n = 50, k = 3, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>());

    ASSERT_EQ(result.info, 0);
    for (int i = 0; i < k; ++i) {
        Real expected = static_cast<Real>(n - i);  // 50, 49, 48
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v))
            << "sigma[" << i << "]";
    }
}

// ============================================================================
// lansvd: wide matrix (m < n)
// ============================================================================

TYPED_TEST(LansvdTest, WideMatrix) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 80, n = 120, k = 3, kmax = 50;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>());

    ASSERT_EQ(result.info, 0);
    for (int i = 0; i < k; ++i) {
        Real expected = static_cast<Real>(m - i);  // 80, 79, 78
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v))
            << "sigma[" << i << "]";
    }
}

// ============================================================================
// lansvd: single singular value (k=1)
// ============================================================================

TYPED_TEST(LansvdTest, SingleSingularValue) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 60, n = 40, k = 1, kmax = 20;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>());

    ASSERT_EQ(result.info, 0);
    ASSERT_EQ(static_cast<int>(result.sigma.size()), 1);
    Real expected = static_cast<Real>(n);  // largest = 40
    EXPECT_NEAR(static_cast<double>(result.sigma[0]),
                static_cast<double>(expected),
                static_cast<double>(Tol<Scalar>::v));
}

// ============================================================================
// lansvd: lambda captures external data
// ============================================================================

TYPED_TEST(LansvdTest, LambdaCapture) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 50, n = 50, k = 3, kmax = 50;

    // Store diagonal entries in a captured vector.
    std::vector<Real> diag(n);
    for (int i = 0; i < n; ++i) diag[i] = static_cast<Real>(10 * (i + 1));

    auto result = propack::lansvd<Scalar>(m, n, k, kmax,
        [&diag](char trans, int m, int n, const Scalar* x, Scalar* y) {
            int len = (trans == 'N' || trans == 'n') ? m : n;
            int mn  = std::min(m, n);
            for (int i = 0; i < len; ++i) y[i] = Scalar(0);
            for (int i = 0; i < mn; ++i)
                y[i] = diag[i] * x[i];
        });

    ASSERT_EQ(result.info, 0);
    // Largest singular values: 500, 490, 480
    for (int i = 0; i < k; ++i) {
        Real expected = static_cast<Real>(10 * (n - i));
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v) * 10.0)
            << "sigma[" << i << "]";
    }
}

// ============================================================================
// lansvd: custom options (MGS + extended local reorth)
// ============================================================================

TYPED_TEST(LansvdTest, CustomOptions) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 80, n = 60, k = 3, kmax = 50;

    propack::Options opts;
    opts.use_mgs = true;
    opts.extended_local_reorth = true;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>(),
                                          true, true, 0.0, opts);

    ASSERT_EQ(result.info, 0);
    for (int i = 0; i < k; ++i) {
        Real expected = static_cast<Real>(n - i);
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v));
    }
}

// ============================================================================
// lansvd: A*V(:,i) ≈ sigma_i * U(:,i) residual check
// ============================================================================

TYPED_TEST(LansvdTest, Residual) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, k = 5, kmax = 50;

    auto aprod = diagonal_aprod<Scalar>();
    auto result = propack::lansvd<Scalar>(m, n, k, kmax, aprod, true, true);

    ASSERT_EQ(result.info, 0);

    std::vector<Scalar> Av(m);
    for (int j = 0; j < k; ++j) {
        // Compute A * V(:,j)
        const Scalar* vj = result.V.data() + static_cast<size_t>(j) * n;
        aprod('N', m, n, vj, Av.data());

        // Compare with sigma[j] * U(:,j)
        const Scalar* uj = result.U.data() + static_cast<size_t>(j) * m;
        double max_err = 0.0;
        for (int i = 0; i < m; ++i) {
            auto diff = Av[i] - static_cast<Real>(result.sigma[j]) * uj[i];
            max_err = std::max(max_err, static_cast<double>(std::abs(diff)));
        }
        EXPECT_LT(max_err, static_cast<double>(Tol<Scalar>::v))
            << "Residual for triplet " << j;
    }
}

// ============================================================================
// Typed test fixture for lansvd_irl
// ============================================================================

template<typename Scalar>
class LansvdIrlTest : public ::testing::Test {};

TYPED_TEST_SUITE(LansvdIrlTest, ScalarTypes);

// ============================================================================
// lansvd_irl: largest singular values
// ============================================================================

TYPED_TEST(LansvdIrlTest, LargestSingularValues) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, nwanted = 5, dim = 30, p = 15, maxiter = 100;

    auto result = propack::lansvd_irl<Scalar>('L', m, n, nwanted, dim, p, maxiter,
                                              diagonal_aprod<Scalar>());

    EXPECT_EQ(result.info, 0);
    EXPECT_EQ(static_cast<int>(result.sigma.size()), nwanted);

    for (int i = 0; i < nwanted; ++i) {
        Real expected = static_cast<Real>(n - i);  // 80, 79, ...
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v))
            << "sigma[" << i << "]";
    }
}

// ============================================================================
// lansvd_irl: smallest singular values
// ============================================================================

TYPED_TEST(LansvdIrlTest, SmallestSingularValues) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, nwanted = 3, dim = 30, p = 15, maxiter = 200;

    auto result = propack::lansvd_irl<Scalar>('S', m, n, nwanted, dim, p, maxiter,
                                              diagonal_aprod<Scalar>());

    EXPECT_EQ(result.info, 0);
    EXPECT_EQ(static_cast<int>(result.sigma.size()), nwanted);

    // Smallest singular values of diag(1,...,80): 1, 2, 3
    for (int i = 0; i < nwanted; ++i) {
        Real expected = static_cast<Real>(nwanted - i);  // 3, 2, 1
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v))
            << "sigma[" << i << "]";
    }
}

// ============================================================================
// lansvd_irl: no vectors mode
// ============================================================================

TYPED_TEST(LansvdIrlTest, NoVectors) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, nwanted = 5, dim = 30, p = 15, maxiter = 100;

    auto result = propack::lansvd_irl<Scalar>('L', m, n, nwanted, dim, p, maxiter,
                                              diagonal_aprod<Scalar>(),
                                              false, false);

    EXPECT_EQ(result.info, 0);
    EXPECT_TRUE(result.U.empty());
    EXPECT_TRUE(result.V.empty());
    for (int i = 0; i < nwanted; ++i) {
        Real expected = static_cast<Real>(n - i);
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v));
    }
}

// ============================================================================
// lansvd_irl: orthonormality of U and V
// ============================================================================

TYPED_TEST(LansvdIrlTest, VectorsOrthonormal) {
    using Scalar = TypeParam;

    const int m = 100, n = 80, nwanted = 5, dim = 30, p = 15, maxiter = 100;

    auto result = propack::lansvd_irl<Scalar>('L', m, n, nwanted, dim, p, maxiter,
                                              diagonal_aprod<Scalar>(),
                                              true, true);

    ASSERT_EQ(result.info, 0);
    ASSERT_EQ(static_cast<int>(result.U.size()), m * nwanted);
    ASSERT_EQ(static_cast<int>(result.V.size()), n * nwanted);

    double u_err = orthogonality_error(result.U.data(), m, nwanted);
    double v_err = orthogonality_error(result.V.data(), n, nwanted);
    EXPECT_LT(u_err, static_cast<double>(Tol<Scalar>::v)) << "U not orthonormal";
    EXPECT_LT(v_err, static_cast<double>(Tol<Scalar>::v)) << "V not orthonormal";
}

// ============================================================================
// lansvd_irl: residual check A*V(:,i) ≈ sigma_i * U(:,i)
// ============================================================================

TYPED_TEST(LansvdIrlTest, Residual) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, nwanted = 5, dim = 30, p = 15, maxiter = 100;

    auto aprod = diagonal_aprod<Scalar>();
    auto result = propack::lansvd_irl<Scalar>('L', m, n, nwanted, dim, p, maxiter,
                                              aprod, true, true);

    ASSERT_EQ(result.info, 0);

    std::vector<Scalar> Av(m);
    for (int j = 0; j < nwanted; ++j) {
        const Scalar* vj = result.V.data() + static_cast<size_t>(j) * n;
        aprod('N', m, n, vj, Av.data());

        const Scalar* uj = result.U.data() + static_cast<size_t>(j) * m;
        double max_err = 0.0;
        for (int i = 0; i < m; ++i) {
            auto diff = Av[i] - static_cast<Real>(result.sigma[j]) * uj[i];
            max_err = std::max(max_err, static_cast<double>(std::abs(diff)));
        }
        EXPECT_LT(max_err, static_cast<double>(Tol<Scalar>::v))
            << "Residual for triplet " << j;
    }
}

// ============================================================================
// lansvd_irl: custom IRL options (min_relgap)
// ============================================================================

TYPED_TEST(LansvdIrlTest, CustomOptions) {
    using Scalar = TypeParam;
    using Real   = typename RealOf<Scalar>::type;

    const int m = 100, n = 80, nwanted = 3, dim = 25, p = 12, maxiter = 100;

    propack::IRLOptions opts;
    opts.use_mgs = true;
    opts.min_relgap = 1e-6;

    auto result = propack::lansvd_irl<Scalar>('L', m, n, nwanted, dim, p, maxiter,
                                              diagonal_aprod<Scalar>(),
                                              true, true, 0.0, opts);

    ASSERT_EQ(result.info, 0);
    for (int i = 0; i < nwanted; ++i) {
        Real expected = static_cast<Real>(n - i);
        EXPECT_NEAR(static_cast<double>(result.sigma[i]),
                    static_cast<double>(expected),
                    static_cast<double>(Tol<Scalar>::v));
    }
}

// ============================================================================
// SVDResult metadata
// ============================================================================

TYPED_TEST(LansvdTest, ResultMetadata) {
    using Scalar = TypeParam;

    const int m = 80, n = 60, k = 4, kmax = 30;

    auto result = propack::lansvd<Scalar>(m, n, k, kmax, diagonal_aprod<Scalar>());

    EXPECT_EQ(result.m, m);
    EXPECT_EQ(result.n, n);
    EXPECT_EQ(result.k, k);
}

TYPED_TEST(LansvdIrlTest, ResultMetadata) {
    using Scalar = TypeParam;

    const int m = 80, n = 60, nwanted = 4, dim = 25, p = 12, maxiter = 50;

    auto result = propack::lansvd_irl<Scalar>('L', m, n, nwanted, dim, p, maxiter,
                                              diagonal_aprod<Scalar>());

    EXPECT_EQ(result.m, m);
    EXPECT_EQ(result.n, n);
    EXPECT_EQ(result.k, nwanted);
}
