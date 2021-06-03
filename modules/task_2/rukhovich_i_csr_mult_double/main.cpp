// Copyright 2021 Igor Rukhovich
#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "../../modules/task_2/rukhovich_i_csr_mult_double/csr_mult_double.h"

template <class Matrix>
double MeasureMulTime(const Matrix lhs, const Matrix rhs, size_t num_threads = 1,
                      size_t num_loops = 15) {
  omp_set_num_threads(num_threads);
  std::vector<double> times;
  times.reserve(num_loops);
  double ans = 0;
  for (size_t i = 0; i < num_loops; ++i) {
    Matrix lhs_tmp = lhs, rhs_tmp = rhs;
    double start = omp_get_wtime();
    lhs_tmp *= rhs_tmp;
    times.emplace_back(omp_get_wtime() - start);
    ans += times.back();
  }
  std::sort(times.begin(), times.end());
  if (times.size() > 2) {
    ans -= times.front() + times.back();
    return ans / static_cast<double>(times.size() - 2);
  }
  return ans / static_cast<double>(times.size());
}

TEST(OMP, TestMultEmpty) {
  CSRMatrixOMP<double> l_dmat(2, 5), r_dmat(5, 10);
  l_dmat *= r_dmat;
  CSRMatrixOMP<double> res_dmat(2, 10);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(OMP, TestMultDiag) {
  std::vector<double> diag_mat(10 * 10);
  for (size_t i = 0; i < 10; ++i) {
    diag_mat[i * 10 + i] = 1;
  }
  CSRMatrixOMP<double> l_dmat(10, 10, diag_mat), r_dmat(10, 10, diag_mat);
  l_dmat *= r_dmat;
  CSRMatrixOMP<double> res_dmat(10, 10, diag_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(OMP, TestMultSimpleToDiag) {
  std::vector<double> diag_mat(10 * 10), init_mat(10 * 10);
  for (size_t i = 0; i < 10; ++i) {
    diag_mat[i * 10 + i] = 1;
    for (size_t j = 0; j < 10; ++j) {
      init_mat[i * 10 + j] = RandomDouble::Next();
    }
  }
  CSRMatrixOMP<double> l_dmat(10, 10, init_mat), r_dmat(10, 10, diag_mat);
  l_dmat *= r_dmat;
  CSRMatrixOMP<double> res_dmat(10, 10, init_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(OMP, TestMultDense) {
  size_t same = 64;
  size_t a = same, b = same, c = same;
  std::vector<double> l_init_mat(a * b), r_init_mat(b * c);
  for (size_t i = 0; i < a; ++i) {
    for (size_t j = 0; j < b; ++j) {
      l_init_mat[i * b + j] = RandomDouble::Next();
    }
  }
  for (size_t j = 0; j < b; ++j) {
    for (size_t k = 0; k < c; ++k) {
      r_init_mat[j * c + k] = RandomDouble::Next();
    }
  }
  std::vector<double> res_mat(a * c);
  for (size_t i = 0; i < a; ++i)
    for (size_t j = 0; j < b; ++j)
      for (size_t k = 0; k < c; ++k) {
        res_mat[i * c + k] += l_init_mat[i * b + j] * r_init_mat[j * c + k];
      }
  CSRMatrixOMP<double> l_dmat(a, b, l_init_mat), r_dmat(b, c, r_init_mat);
  double t1 = MeasureMulTime(l_dmat, r_dmat, 1);
  double t2 = MeasureMulTime(l_dmat, r_dmat, 2);
  double t4 = MeasureMulTime(l_dmat, r_dmat, 4);
  double t8 = MeasureMulTime(l_dmat, r_dmat, 8);
  double t16 = MeasureMulTime(l_dmat, r_dmat, 16);
  std::cout.precision(3);
  std::cout << "Threads:\t1\t2\t4\t8\t16\n";
  std::cout << std::fixed << "Times:\t\t" << t1 << '\t' << t2 << '\t' << t4 << '\t' << t8 << '\t'
            << t16 << '\n';
  std::cout << std::fixed << "Speed up:\t" << 1. << '\t' << t1 / t2 << '\t' << t1 / t4 << '\t'
            << t1 / t8 << '\t' << t1 / t16 << '\n';
  l_dmat *= r_dmat;
  CSRMatrixOMP<double> res_dmat(a, c, res_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(OMP, TestMultSparse) {
  size_t same = 64;
  double zero_per_elem = 2;
  size_t a = same, b = same, c = same;
  std::random_device device;
  std::mt19937 gen(device());
  std::bernoulli_distribution dist(1. / (zero_per_elem + 1.));
  std::vector<double> l_init_mat(a * b), r_init_mat(b * c);
  for (size_t i = 0; i < a; ++i) {
    for (size_t j = 0; j < b; ++j) {
      l_init_mat[i * b + j] = (dist(gen) ? RandomDouble::Next() : 0);
    }
  }
  for (size_t j = 0; j < b; ++j) {
    for (size_t k = 0; k < c; ++k) {
      r_init_mat[j * c + k] = (dist(gen) ? RandomDouble::Next() : 0);
    }
  }
  std::vector<double> res_mat(a * c);
  for (size_t i = 0; i < a; ++i)
    for (size_t j = 0; j < b; ++j)
      for (size_t k = 0; k < c; ++k) {
        res_mat[i * c + k] += l_init_mat[i * b + j] * r_init_mat[j * c + k];
      }
  CSRMatrixOMP<double> l_dmat(a, b, l_init_mat), r_dmat(b, c, r_init_mat);
  double t1 = MeasureMulTime(l_dmat, r_dmat, 1);
  double t2 = MeasureMulTime(l_dmat, r_dmat, 2);
  double t4 = MeasureMulTime(l_dmat, r_dmat, 4);
  double t8 = MeasureMulTime(l_dmat, r_dmat, 8);
  double t16 = MeasureMulTime(l_dmat, r_dmat, 16);
  std::cout.precision(3);
  std::cout << "Threads:\t1\t2\t4\t8\t16\n";
  std::cout << std::fixed << "Times:\t\t" << t1 << '\t' << t2 << '\t' << t4 << '\t' << t8 << '\t'
            << t16 << '\n';
  std::cout << std::fixed << "Speed up:\t" << 1. << '\t' << t1 / t2 << '\t' << t1 / t4 << '\t'
            << t1 / t8 << '\t' << t1 / t16 << '\n';
  l_dmat *= r_dmat;
  CSRMatrixOMP<double> res_dmat(a, c, res_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}
