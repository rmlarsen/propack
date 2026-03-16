// propack.hpp — Header-only C++17 template API for PROPACK SVD drivers
//
// Usage:
//   #include <propack.hpp>
//   auto result = propack::lansvd<double>(m, n, k, kmax,
//       [&](char trans, int m, int n, const double* x, double* y) {
//           // matrix-vector product
//       });
//
// Link against the PROPACK Fortran libraries (e.g. dpropack, spropack, etc.)

#ifndef PROPACK_HPP_
#define PROPACK_HPP_

#include <algorithm>
#include <complex>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace propack {

// ============================================================================
// Options
// ============================================================================

struct Options {
    double delta = 0;              // level of orthogonality (0 = sqrt(eps))
    double eta = 0;                // reorth threshold (0 = eps^(3/4)*3/4)
    double anorm = 0;              // ||A|| estimate (0 = auto)
    double eps1_aprod = 0;         // rounding error for A*v (0 = auto: 5*sqrt(n)*eps)
    double eps1_atprod = 0;        // rounding error for A'*u (0 = auto: 5*sqrt(m)*eps)
    bool use_mgs = false;          // use iterated MGS instead of CGS
    bool extended_local_reorth = false;  // extended local reorthogonalization
};

struct IRLOptions : Options {
    double min_relgap = 1e-3;      // smallest relative gap between shifts
};

// ============================================================================
// SVDResult
// ============================================================================

namespace detail {
    template<typename T> struct real_type                   { using type = T; };
    template<typename T> struct real_type<std::complex<T>>  { using type = T; };
    template<typename T> using real_type_t = typename real_type<T>::type;
}

template<typename Scalar>
struct SVDResult {
    using Real = detail::real_type_t<Scalar>;

    std::vector<Real> sigma;   // singular values
    std::vector<Real> bnd;     // error bounds
    std::vector<Scalar> U;     // left singular vectors (column-major)
    std::vector<Scalar> V;     // right singular vectors (column-major)
    int m = 0, n = 0, k = 0;
    int info = 0;              // 0 = success, >0 invariant subspace, <0 partial convergence
};

// ============================================================================
// APROD callback type
// ============================================================================

template<typename Scalar>
using AprodFn = std::function<void(char, int, int, const Scalar*, Scalar*)>;

// ============================================================================
// extern "C" declarations for Fortran bind(C) wrappers
// ============================================================================

namespace detail {

// C callback type: void(char, int, int, const Scalar*, Scalar*, void*)
template<typename Scalar>
using c_aprod_fn = void(*)(char, int, int, const Scalar*, Scalar*, void*);

extern "C" {

// Double precision
void propack_dlansvd(char jobu, char jobv, int m, int n, int k, int kmax,
                     c_aprod_fn<double> aprod, void* userdata,
                     double* U, int ldu, double* Sigma, double* bnd,
                     double* V, int ldv, double tolin,
                     const double* doption, const int* ioption, int* info);
void propack_dlansvd_irl(char which, char jobu, char jobv,
                         int m, int n, int dim, int p, int nwanted,
                         int maxiter,
                         c_aprod_fn<double> aprod, void* userdata,
                         double* U, int ldu, double* Sigma, double* bnd,
                         double* V, int ldv, double tolin,
                         const double* doption, const int* ioption,
                         int* info);

// Single precision
void propack_slansvd(char jobu, char jobv, int m, int n, int k, int kmax,
                     c_aprod_fn<float> aprod, void* userdata,
                     float* U, int ldu, float* Sigma, float* bnd,
                     float* V, int ldv, float tolin,
                     const float* doption, const int* ioption, int* info);
void propack_slansvd_irl(char which, char jobu, char jobv,
                         int m, int n, int dim, int p, int nwanted,
                         int maxiter,
                         c_aprod_fn<float> aprod, void* userdata,
                         float* U, int ldu, float* Sigma, float* bnd,
                         float* V, int ldv, float tolin,
                         const float* doption, const int* ioption,
                         int* info);

// Complex double precision
void propack_zlansvd(char jobu, char jobv, int m, int n, int k, int kmax,
                     c_aprod_fn<std::complex<double>> aprod, void* userdata,
                     std::complex<double>* U, int ldu,
                     double* Sigma, double* bnd,
                     std::complex<double>* V, int ldv, double tolin,
                     const double* doption, const int* ioption, int* info);
void propack_zlansvd_irl(char which, char jobu, char jobv,
                         int m, int n, int dim, int p, int nwanted,
                         int maxiter,
                         c_aprod_fn<std::complex<double>> aprod,
                         void* userdata,
                         std::complex<double>* U, int ldu,
                         double* Sigma, double* bnd,
                         std::complex<double>* V, int ldv, double tolin,
                         const double* doption, const int* ioption,
                         int* info);

// Complex single precision
void propack_clansvd(char jobu, char jobv, int m, int n, int k, int kmax,
                     c_aprod_fn<std::complex<float>> aprod, void* userdata,
                     std::complex<float>* U, int ldu,
                     float* Sigma, float* bnd,
                     std::complex<float>* V, int ldv, float tolin,
                     const float* doption, const int* ioption, int* info);
void propack_clansvd_irl(char which, char jobu, char jobv,
                         int m, int n, int dim, int p, int nwanted,
                         int maxiter,
                         c_aprod_fn<std::complex<float>> aprod,
                         void* userdata,
                         std::complex<float>* U, int ldu,
                         float* Sigma, float* bnd,
                         std::complex<float>* V, int ldv, float tolin,
                         const float* doption, const int* ioption,
                         int* info);

}  // extern "C"

// ============================================================================
// Trampoline: cast void* back to std::function and invoke
// ============================================================================

template<typename Scalar>
inline void trampoline(char t, int m, int n,
                       const Scalar* x, Scalar* y, void* ud) {
    (*static_cast<AprodFn<Scalar>*>(ud))(t, m, n, x, y);
}

// ============================================================================
// call_lansvd — dispatch to correct Fortran wrapper via if constexpr
// ============================================================================

template<typename Scalar>
void call_lansvd(char jobu, char jobv, int m, int n, int k, int kmax,
                 AprodFn<Scalar>& fn,
                 Scalar* U, int ldu, real_type_t<Scalar>* sigma,
                 real_type_t<Scalar>* bnd,
                 Scalar* V, int ldv, real_type_t<Scalar> tolin,
                 const real_type_t<Scalar>* doption, const int* ioption,
                 int* info) {
    using S = Scalar;
    if constexpr (std::is_same_v<S, double>)
        propack_dlansvd(jobu, jobv, m, n, k, kmax,
                        trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                        tolin, doption, ioption, info);
    else if constexpr (std::is_same_v<S, float>)
        propack_slansvd(jobu, jobv, m, n, k, kmax,
                        trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                        tolin, doption, ioption, info);
    else if constexpr (std::is_same_v<S, std::complex<double>>)
        propack_zlansvd(jobu, jobv, m, n, k, kmax,
                        trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                        tolin, doption, ioption, info);
    else if constexpr (std::is_same_v<S, std::complex<float>>)
        propack_clansvd(jobu, jobv, m, n, k, kmax,
                        trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                        tolin, doption, ioption, info);
}

// ============================================================================
// call_lansvd_irl — dispatch to correct Fortran wrapper via if constexpr
// ============================================================================

template<typename Scalar>
void call_lansvd_irl(char which, char jobu, char jobv,
                     int m, int n, int dim, int p, int nwanted, int maxiter,
                     AprodFn<Scalar>& fn,
                     Scalar* U, int ldu, real_type_t<Scalar>* sigma,
                     real_type_t<Scalar>* bnd,
                     Scalar* V, int ldv, real_type_t<Scalar> tolin,
                     const real_type_t<Scalar>* doption, const int* ioption,
                     int* info) {
    using S = Scalar;
    if constexpr (std::is_same_v<S, double>)
        propack_dlansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter,
                            trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                            tolin, doption, ioption, info);
    else if constexpr (std::is_same_v<S, float>)
        propack_slansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter,
                            trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                            tolin, doption, ioption, info);
    else if constexpr (std::is_same_v<S, std::complex<double>>)
        propack_zlansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter,
                            trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                            tolin, doption, ioption, info);
    else if constexpr (std::is_same_v<S, std::complex<float>>)
        propack_clansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter,
                            trampoline<S>, &fn, U, ldu, sigma, bnd, V, ldv,
                            tolin, doption, ioption, info);
}

// ============================================================================
// Option conversion helpers
// ============================================================================

template<typename Real>
inline void fill_doption(Real* doption, const Options& opts) {
    doption[0] = static_cast<Real>(opts.delta);
    doption[1] = static_cast<Real>(opts.eta);
    doption[2] = static_cast<Real>(opts.anorm);
    doption[3] = Real(0);  // unused by lansvd (only by lansvd_irl)
    doption[4] = static_cast<Real>(opts.eps1_aprod);
    doption[5] = static_cast<Real>(opts.eps1_atprod);
}

template<typename Real>
inline void fill_doption_irl(Real* doption, const IRLOptions& opts) {
    fill_doption(doption, opts);
    doption[3] = static_cast<Real>(opts.min_relgap);
    // doption[4] and doption[5] already set by fill_doption
}

inline void fill_ioption(int* ioption, const Options& opts) {
    ioption[0] = opts.use_mgs ? 1 : 0;
    ioption[1] = opts.extended_local_reorth ? 1 : 0;
}

}  // namespace detail

// ============================================================================
// propack::lansvd — Lanczos bidiagonalization SVD
// ============================================================================

template<typename Scalar>
[[nodiscard]] SVDResult<Scalar> lansvd(
        int m, int n, int k, int kmax,
        AprodFn<Scalar> aprod,
        bool compute_u = true, bool compute_v = true,
        double tol = 0.0,
        const Options& opts = {}) {
    using Real = detail::real_type_t<Scalar>;

    char jobu = compute_u ? 'Y' : 'N';
    char jobv = compute_v ? 'Y' : 'N';

    // Allocate Fortran-sized buffers (kmax+1 columns for U, kmax for V)
    std::vector<Scalar> U_buf(static_cast<size_t>(m) * (kmax + 1), Scalar(0));
    std::vector<Scalar> V_buf(static_cast<size_t>(n) * kmax, Scalar(0));
    std::vector<Real> sigma(k);
    std::vector<Real> bnd(k);

    Real doption[6];
    int ioption[2];
    detail::fill_doption(doption, opts);
    detail::fill_ioption(ioption, opts);

    int info = 0;
    detail::call_lansvd(jobu, jobv, m, n, k, kmax, aprod,
                        U_buf.data(), m, sigma.data(), bnd.data(),
                        V_buf.data(), n,
                        static_cast<Real>(tol),
                        doption, ioption, &info);

    SVDResult<Scalar> result;
    result.m = m;
    result.n = n;
    result.k = k;
    result.info = info;
    result.sigma = std::move(sigma);
    result.bnd = std::move(bnd);

    if (compute_u) {
        result.U.resize(static_cast<size_t>(m) * k);
        for (int j = 0; j < k; ++j)
            std::copy(U_buf.data() + static_cast<size_t>(j) * m,
                      U_buf.data() + static_cast<size_t>(j) * m + m,
                      result.U.data() + static_cast<size_t>(j) * m);
    }
    if (compute_v) {
        result.V.resize(static_cast<size_t>(n) * k);
        for (int j = 0; j < k; ++j)
            std::copy(V_buf.data() + static_cast<size_t>(j) * n,
                      V_buf.data() + static_cast<size_t>(j) * n + n,
                      result.V.data() + static_cast<size_t>(j) * n);
    }

    return result;
}

// ============================================================================
// propack::lansvd_irl — Implicitly restarted Lanczos SVD
// ============================================================================

template<typename Scalar>
[[nodiscard]] SVDResult<Scalar> lansvd_irl(
        char which, int m, int n,
        int nwanted, int dim, int p, int maxiter,
        AprodFn<Scalar> aprod,
        bool compute_u = true, bool compute_v = true,
        double tol = 0.0,
        const IRLOptions& opts = {}) {
    using Real = detail::real_type_t<Scalar>;

    char jobu = compute_u ? 'Y' : 'N';
    char jobv = compute_v ? 'Y' : 'N';

    // Allocate Fortran-sized buffers
    std::vector<Scalar> U_buf(static_cast<size_t>(m) * (dim + 1), Scalar(0));
    std::vector<Scalar> V_buf(static_cast<size_t>(n) * dim, Scalar(0));
    std::vector<Real> sigma(nwanted);
    std::vector<Real> bnd(nwanted);

    Real doption[6];
    int ioption[2];
    detail::fill_doption_irl(doption, opts);
    detail::fill_ioption(ioption, opts);

    int info = 0;
    detail::call_lansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter,
                            aprod, U_buf.data(), m, sigma.data(), bnd.data(),
                            V_buf.data(), n,
                            static_cast<Real>(tol),
                            doption, ioption, &info);

    SVDResult<Scalar> result;
    result.m = m;
    result.n = n;
    result.k = nwanted;
    result.info = info;
    result.sigma = std::move(sigma);
    result.bnd = std::move(bnd);

    if (compute_u) {
        result.U.resize(static_cast<size_t>(m) * nwanted);
        for (int j = 0; j < nwanted; ++j)
            std::copy(U_buf.data() + static_cast<size_t>(j) * m,
                      U_buf.data() + static_cast<size_t>(j) * m + m,
                      result.U.data() + static_cast<size_t>(j) * m);
    }
    if (compute_v) {
        result.V.resize(static_cast<size_t>(n) * nwanted);
        for (int j = 0; j < nwanted; ++j)
            std::copy(V_buf.data() + static_cast<size_t>(j) * n,
                      V_buf.data() + static_cast<size_t>(j) * n + n,
                      result.V.data() + static_cast<size_t>(j) * n);
    }

    return result;
}

}  // namespace propack

#endif  // PROPACK_HPP_
