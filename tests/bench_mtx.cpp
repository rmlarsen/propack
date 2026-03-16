// bench_fortran_propack.cpp — Benchmark Fortran PROPACK on Matrix Market
// matrices using the C++ API wrapper.
//
// Usage:  bench_fortran_propack <matrix_dir> [k1 k2 ...]
//         Default k values: 10 50
//
// Links against the Fortran PROPACK dpropack library.
//
// (C) Rasmus Munk Larsen, 2026

// Include Fortran PROPACK C++ wrapper by absolute path to avoid
// collision with propack++ propack.hpp
#include "propack.hpp"

#include <mmio.hpp>

#include <Eigen/Sparse>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct BenchResult {
  std::string matrix;
  std::string algorithm;
  int k;
  Eigen::Index m, n;
  long long nnz;
  int nconv;
  double time_s;
  double max_sigma;
  double min_sigma;
};

void print_json(const std::vector<BenchResult>& results) {
  std::cout << "[\n";
  for (size_t i = 0; i < results.size(); ++i) {
    const auto& r = results[i];
    std::cout << "  {"
              << "\"matrix\":\"" << r.matrix << "\","
              << "\"algorithm\":\"" << r.algorithm << "\","
              << "\"k\":" << r.k << ","
              << "\"m\":" << r.m << ","
              << "\"n\":" << r.n << ","
              << "\"nnz\":" << r.nnz << ","
              << "\"nconv\":" << r.nconv << ","
              << "\"time_s\":" << std::setprecision(6) << r.time_s << ","
              << "\"max_sigma\":" << std::setprecision(12) << r.max_sigma << ","
              << "\"min_sigma\":" << std::setprecision(12) << r.min_sigma
              << "}";
    if (i + 1 < results.size()) std::cout << ",";
    std::cout << "\n";
  }
  std::cout << "]\n";
}

// Sparse matrix-vector product callback for Fortran PROPACK
struct SparseAprod {
  const Eigen::SparseMatrix<double>& A;

  void operator()(char trans, int /*m*/, int /*n*/, const double* x,
                  double* y) const {
    if (trans == 'n' || trans == 'N') {
      Eigen::Map<const Eigen::VectorXd> xv(x, A.cols());
      Eigen::Map<Eigen::VectorXd> yv(y, A.rows());
      yv.noalias() = A * xv;
    } else {
      Eigen::Map<const Eigen::VectorXd> xv(x, A.rows());
      Eigen::Map<Eigen::VectorXd> yv(y, A.cols());
      yv.noalias() = A.adjoint() * xv;
    }
  }
};

BenchResult run_lansvd(const Eigen::SparseMatrix<double>& A,
                       const std::string& name, int k, bool use_mgs) {
  int m = static_cast<int>(A.rows());
  int n = static_cast<int>(A.cols());
  int kmax = std::min(std::min(m, n) - 1, std::max(10 * k + 100, k + 500));

  SparseAprod aprod_obj{A};
  propack::AprodFn<double> aprod = [&](char t, int mm, int nn, const double* x,
                                       double* y) {
    aprod_obj(t, mm, nn, x, y);
  };

  propack::Options opts;
  opts.use_mgs = use_mgs;

  // Warm up
  auto warmup = propack::lansvd<double>(m, n, std::min(k, 5),
                                        std::min(kmax, 50), aprod,
                                        true, true, 0.0, opts);

  // Timed run
  int nruns = (m < 10000 && n < 10000) ? 3 : 1;
  double total_time = 0;
  propack::SVDResult<double> result;
  for (int run = 0; run < nruns; ++run) {
    auto t0 = std::chrono::high_resolution_clock::now();
    result = propack::lansvd<double>(m, n, k, kmax, aprod,
                                     true, true, 0.0, opts);
    auto t1 = std::chrono::high_resolution_clock::now();
    total_time += std::chrono::duration<double>(t1 - t0).count();
  }

  int nconv = (result.info >= 0) ? k : -result.info;
  if (result.info > 0) nconv = k;  // invariant subspace found

  std::string gs_tag = use_mgs ? "_mgs" : "_cgs";
  BenchResult br;
  br.matrix = name;
  br.algorithm = "fortran_lansvd" + gs_tag;
  br.k = k;
  br.m = m;
  br.n = n;
  br.nnz = A.nonZeros();
  br.nconv = nconv;
  br.time_s = total_time / nruns;
  br.max_sigma =
      result.sigma.size() > 0 ? result.sigma[0] : 0;
  br.min_sigma =
      result.sigma.size() > 0 ? result.sigma[result.sigma.size() - 1] : 0;
  return br;
}

BenchResult run_lansvd_irl(const Eigen::SparseMatrix<double>& A,
                           const std::string& name, int k, bool use_mgs) {
  int m = static_cast<int>(A.rows());
  int n = static_cast<int>(A.cols());
  int dim = std::min(std::min(m, n) - 1, std::max(3 * k + 30, k + 50));
  int p = dim - k - 2;
  if (p < 2) {
    p = 2;
    dim = k + p + 2;
  }

  SparseAprod aprod_obj{A};
  propack::AprodFn<double> aprod = [&](char t, int mm, int nn, const double* x,
                                       double* y) {
    aprod_obj(t, mm, nn, x, y);
  };

  propack::IRLOptions opts;
  opts.use_mgs = use_mgs;

  // Warm up
  auto warmup = propack::lansvd_irl<double>(
      'L', m, n, std::min(k, 5), std::min(dim, std::min(k, 5) + 20),
      std::min(p, 10), 50, aprod, true, true, 0.0, opts);

  // Timed run
  int nruns = (m < 10000 && n < 10000) ? 3 : 1;
  double total_time = 0;
  propack::SVDResult<double> result;
  for (int run = 0; run < nruns; ++run) {
    auto t0 = std::chrono::high_resolution_clock::now();
    result = propack::lansvd_irl<double>('L', m, n, k, dim, p, 500, aprod,
                                         true, true, 0.0, opts);
    auto t1 = std::chrono::high_resolution_clock::now();
    total_time += std::chrono::duration<double>(t1 - t0).count();
  }

  int nconv = (result.info >= 0) ? k : -result.info;
  if (result.info > 0) nconv = k;

  std::string gs_tag = use_mgs ? "_mgs" : "_cgs";
  BenchResult br;
  br.matrix = name;
  br.algorithm = "fortran_lansvd_irl" + gs_tag;
  br.k = k;
  br.m = m;
  br.n = n;
  br.nnz = A.nonZeros();
  br.nconv = nconv;
  br.time_s = total_time / nruns;
  br.max_sigma =
      result.sigma.size() > 0 ? result.sigma[0] : 0;
  br.min_sigma =
      result.sigma.size() > 0 ? result.sigma[result.sigma.size() - 1] : 0;
  return br;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <matrix_dir> [k1 k2 ...]\n";
    return 1;
  }

  std::string matrix_dir = argv[1];
  std::vector<int> k_values;
  for (int i = 2; i < argc; ++i) {
    k_values.push_back(std::atoi(argv[i]));
  }
  if (k_values.empty()) {
    k_values = {10, 50};
  }

  // Collect all .mtx files
  std::vector<std::pair<std::string, std::string>> matrices;
  for (const auto& entry : fs::directory_iterator(matrix_dir)) {
    if (entry.path().extension() == ".mtx") {
      std::string name = entry.path().stem().string();
      matrices.emplace_back(name, entry.path().string());
    }
  }
  std::sort(matrices.begin(), matrices.end());

  std::vector<BenchResult> all_results;

  for (const auto& [name, path] : matrices) {
    std::cerr << "Loading " << name << " from " << path << "... "
              << std::flush;
    try {
      auto A = propack::read_matrix_market<double>(path);
      std::cerr << A.rows() << "x" << A.cols() << " nnz=" << A.nonZeros()
                << "\n";

      // Skip apache2 — IRL doesn't converge within reasonable time
      if (name == "apache2") {
        std::cerr << "  Skipping (too slow for Fortran PROPACK IRL)\n";
        continue;
      }

      for (int k : k_values) {
        if (k >= std::min(A.rows(), A.cols()) - 1) {
          std::cerr << "  Skipping k=" << k << " (too large for matrix)\n";
          continue;
        }

        for (bool use_mgs : {false, true}) {
          std::string gs = use_mgs ? "MGS" : "CGS";

          std::cerr << "  fortran_lansvd(" << gs << ") k=" << k << "... "
                    << std::flush;
          auto r1 = run_lansvd(A, name, k, use_mgs);
          std::cerr << r1.time_s << "s nconv=" << r1.nconv << "\n";
          all_results.push_back(r1);

          std::cerr << "  fortran_lansvd_irl(" << gs << ") k=" << k << "... "
                    << std::flush;
          auto r2 = run_lansvd_irl(A, name, k, use_mgs);
          std::cerr << r2.time_s << "s nconv=" << r2.nconv << "\n";
          all_results.push_back(r2);
        }
      }
    } catch (const std::exception& e) {
      std::cerr << "ERROR: " << e.what() << "\n";
    }
  }

  print_json(all_results);
  return 0;
}
