// mmio.hpp — Simple Matrix Market reader for Eigen sparse matrices
//
// (C) Rasmus Munk Larsen, 2026

#ifndef PROPACK_MMIO_HPP_
#define PROPACK_MMIO_HPP_

#include <Eigen/SparseCore>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace propack {

/// Read a Matrix Market file into an Eigen sparse matrix.
/// Supports: coordinate real/pattern general/symmetric.
template <typename Scalar = double>
Eigen::SparseMatrix<Scalar> read_matrix_market(const std::string& filename) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  // Parse header line
  std::string line;
  std::getline(in, line);

  bool is_symmetric = (line.find("symmetric") != std::string::npos);
  bool is_pattern = (line.find("pattern") != std::string::npos);
  bool is_complex = (line.find("complex") != std::string::npos);

  if (is_complex) {
    throw std::runtime_error("Complex Matrix Market files not supported");
  }

  // Skip comments
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '%') continue;
    break;
  }

  // Parse dimensions
  std::istringstream dim_stream(line);
  Eigen::Index m, n;
  long long nnz;
  dim_stream >> m >> n >> nnz;

  // Read entries
  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(is_symmetric ? 2 * nnz : nnz);

  for (long long k = 0; k < nnz; ++k) {
    Eigen::Index i, j;
    Scalar val = Scalar(1);
    in >> i >> j;
    if (!is_pattern) {
      double v;
      in >> v;
      val = static_cast<Scalar>(v);
    }
    --i;  // 1-based to 0-based
    --j;
    triplets.emplace_back(i, j, val);
    if (is_symmetric && i != j) {
      triplets.emplace_back(j, i, val);
    }
  }

  Eigen::SparseMatrix<Scalar> A(m, n);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  return A;
}

}  // namespace propack

#endif  // PROPACK_MMIO_HPP_
