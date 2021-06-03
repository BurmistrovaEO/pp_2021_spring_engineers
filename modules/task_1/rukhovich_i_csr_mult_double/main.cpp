// Copyright 2021 Igor Rukhovich
#include <gtest/gtest.h>

#include <vector>

#include "../../modules/task_1/rukhovich_i_csr_mult_double/csr_mult_double.h"

TEST(Sequential, TestCreate) {
  CSRMatrix<int> mat(0, 0);
  CSRMatrix<double> dmat(0, 0);
  CSRMatrix<double> little_dmat(5, 5);
  std::vector<double> empty_vec(5 * 5, 0);
  CSRMatrix<double> little_dmat_with_zeros(5, 5, empty_vec);
  ASSERT_TRUE(little_dmat == little_dmat_with_zeros);
}

TEST(Sequential, TestMultEmpty) {
  CSRMatrix<double> l_dmat(2, 5), r_dmat(5, 10);
  l_dmat *= r_dmat;
  CSRMatrix<double> res_dmat(2, 10);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(Sequential, TestMultDiag) {
  std::vector<double> diag_mat(10 * 10);
  for (size_t i = 0; i < 10; ++i) {
    diag_mat[i * 10 + i] = 1;
  }
  CSRMatrix<double> l_dmat(10, 10, diag_mat), r_dmat(10, 10, diag_mat);
  l_dmat *= r_dmat;
  CSRMatrix<double> res_dmat(10, 10, diag_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(Sequential, TestMultSimpleToDiag) {
  std::vector<double> diag_mat(10 * 10), init_mat(10 * 10);
  for (size_t i = 0; i < 10; ++i) {
    diag_mat[i * 10 + i] = 1;
    for (size_t j = 0; j < 10; ++j) {
      init_mat[i * 10 + j] = RandomDouble::Next();
    }
  }
  CSRMatrix<double> l_dmat(10, 10, init_mat), r_dmat(10, 10, diag_mat);
  l_dmat *= r_dmat;
  CSRMatrix<double> res_dmat(10, 10, init_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}

TEST(Sequential, TestMult) {
  std::vector<double> l_init_mat(10 * 10), r_init_mat(10 * 10);
  for (size_t i = 0; i < 10; ++i) {
    for (size_t j = 0; j < 10; ++j) {
      l_init_mat[i * 10 + j] = RandomDouble::Next();
      r_init_mat[i * 10 + j] = RandomDouble::Next();
    }
  }
  std::vector<double> res_mat(10 * 10);
  for (size_t i = 0; i < 10; ++i)
    for (size_t j = 0; j < 10; ++j)
      for (size_t k = 0; k < 10; ++k) {
        res_mat[i * 10 + j] += l_init_mat[i * 10 + k] * r_init_mat[k * 10 + j];
      }
  CSRMatrix<double> l_dmat(10, 10, l_init_mat), r_dmat(10, 10, r_init_mat);
  l_dmat *= r_dmat;
  CSRMatrix<double> res_dmat(10, 10, res_mat);
  ASSERT_TRUE(l_dmat == res_dmat);
}
