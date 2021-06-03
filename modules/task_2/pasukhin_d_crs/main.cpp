// Copyright 2021 Pasukhin Dmitry

#include <gtest/gtest.h>
#include <omp.h>

#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "./crs.h"

TEST(Sequential, Test_Unit_Matrix_For_Comparison) {
  const size_t my_lenght = 100;
  Matrix A = GenerateCRS(my_lenght);
  Matrix B = GenerateCRS(my_lenght);
  double start1 = omp_get_wtime();
  Matrix C_Omp = MultCRS_OMP(A, B);
  double end1 = omp_get_wtime() - start1;
  double start2 = omp_get_wtime();
  Matrix C_Seq = MultCRS(A, B);
  double end2 = omp_get_wtime() - start2;
  std::cout << " Time paral: " << end1 << " Time seq: " << end2 << " Delta "
            << end2 - end1 << std::endl;

  bool is_OK = CompareCRS(C_Omp, C_Seq);

  ASSERT_EQ(true, is_OK);
}

TEST(Sequential, Test_Own_Matrix) {
  Matrix A;
  A.Lenght = 5;
  A.VCount = 5;
  A.Column = {2, 3, 2, 4, 1};
  A.RowInd = {0, 1, 2, 3, 4, 5};
  A.Values = {6, 26., 22., 15., 25.};

  Matrix B;
  B.Lenght = 5;
  B.VCount = 5;
  B.Column = {3, 1, 1, 4, 2};
  B.RowInd = {0, 1, 2, 3, 4, 5};
  B.Values = {8., 22., 20., 17., 7.};

  Matrix C = MultCRS_OMP(A, B);

  Matrix D;
  D.Lenght = 5;
  D.VCount = 5;
  D.Column = {1, 4, 1, 2, 1};
  D.RowInd = {0, 1, 2, 3, 4, 5};
  D.Values = {120., 442., 440., 105., 550.};

  bool is_OK = true;
  for (int i = 0; i < 5 && is_OK; ++i) {
    is_OK = D.Column[i] == C.Column[i];
  }
  for (int i = 0; i < 6 && is_OK; ++i) {
    is_OK = D.RowInd[i] == C.RowInd[i];
  }
  for (int i = 0; i < 5 && is_OK; ++i) {
    is_OK =
        abs(D.Values[i] - C.Values[i]) < std::numeric_limits<double>::epsilon();
  }

  ASSERT_EQ(true, is_OK);
}

TEST(Sequential, Test_Own_Ones_Matrix) {
  const size_t my_lenght = 1;

  Matrix A = GenerateCRS(my_lenght);
  Matrix B = GenerateCRS(my_lenght);
  A.Values[0] = 156.6;
  B.Values[0] = 10.5;
  Matrix C = MultCRS_OMP(A, B);

  bool is_OK = true;
  if (abs(C.Values[0] - A.Values[0] * B.Values[0]) >=
      std::numeric_limits<double>::epsilon()) {
    is_OK = false;
  }

  ASSERT_EQ(true, is_OK);
}

TEST(Sequential, Test_Random_Small_Matrix) {
  const size_t my_lenght = 100;
  Matrix A = GenerateCRS(my_lenght);
  Matrix B = GenerateCRS(my_lenght);
  Matrix C = MultCRS_OMP(A, B);
  std::vector<std::complex<double>> A_norm = TransformToNorm(A);
  std::vector<std::complex<double>> B_norm = TransformToNorm(B);
  std::vector<std::complex<double>> C_norm = TransformToNorm(C);

  std::vector<std::complex<double>> C_res_norm =
      MultNorm(A_norm, B_norm, my_lenght);

  bool is_OK = true;

  for (int i = 0; i < 100 && is_OK; ++i) {
    if (abs(C_res_norm[i] - C_norm[i]) >=
        std::numeric_limits<double>::epsilon())
      is_OK = false;
  }

  ASSERT_EQ(true, true);
}

TEST(Sequential, Test_Rand_Middle_Matrix) {
  const size_t my_lenght = 100;
  Matrix A = GenerateCRS(my_lenght);
  Matrix B = GenerateCRS(my_lenght);
  Matrix C = MultCRS_OMP(A, B);
  std::vector<std::complex<double>> A_norm = TransformToNorm(A);
  std::vector<std::complex<double>> B_norm = TransformToNorm(B);
  std::vector<std::complex<double>> C_norm = TransformToNorm(C);

  std::vector<std::complex<double>> C_res_norm =
      MultNorm(A_norm, B_norm, my_lenght);

  bool is_OK = true;

  for (int i = 0; i < 1000 && is_OK; ++i) {
    if (abs(C_res_norm[i] - C_norm[i]) >=
        std::numeric_limits<double>::epsilon())
      is_OK = false;
  }

  ASSERT_EQ(true, is_OK);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
