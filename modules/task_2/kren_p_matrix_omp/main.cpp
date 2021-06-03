// Copyright 2021 Kren Polina
#include <gtest/gtest.h>
#include <omp.h>

#include <iostream>
#include <vector>

#include "../../modules/task_2/kren_p_matrix_omp/matrix.h"

TEST(Strassen_OMP, correct_addition_matrix) {
  int MatrixSize = 20;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> matrix3(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);

  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = 2 * i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = MatrixResultSize - 2 * i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix3[i] = MatrixResultSize;
  }
  AdditionMatrix(matrix1, matrix2, &MatrixResult);
  ASSERT_TRUE(CheckEqualMatrix(matrix3, MatrixResult, MatrixSize));
}

TEST(Strassen_OMP, correct_subtraction_matrix) {
  int MatrixSize = 50;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> matrix3(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = 3 * i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = 3 * i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix3[i] = 0;
  }
  SubtractionMatrix(matrix1, matrix2, &MatrixResult);
  ASSERT_TRUE(CheckEqualMatrix(matrix3, MatrixResult, MatrixSize));
}

TEST(Strassen_OMP, correct_strassen_omp_16x16) {
  int MatrixSize = 16;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> matrix3(MatrixResultSize);
  std::vector<double> matrix4(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = i;
  }
  MultiplicationMatrix(matrix1, matrix2, &matrix3, MatrixSize);
  StrassenAlgorithm_omp(matrix1, matrix2, &matrix4);
  ASSERT_TRUE(CheckEqualMatrix(matrix3, matrix4, MatrixSize));
}

TEST(Strassen_OMP, correct_strassen_omp_64x64) {
  int MatrixSize = 64;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> matrix3(MatrixResultSize);
  std::vector<double> matrix4(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = i;
  }
  MultiplicationMatrix(matrix1, matrix2, &matrix3, MatrixSize);
  StrassenAlgorithm_omp(matrix1, matrix2, &matrix4);
  ASSERT_TRUE(CheckEqualMatrix(matrix3, matrix4, MatrixSize));
}

TEST(Strassen_OMP, correct_strassen_omp_128x128) {
  int MatrixSize = 128;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> matrix3(MatrixResultSize);
  std::vector<double> matrix4(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = MatrixResultSize + i;
  }
  MultiplicationMatrix(matrix1, matrix2, &matrix3, MatrixSize);
  StrassenAlgorithm_omp(matrix1, matrix2, &matrix4);
  ASSERT_TRUE(CheckEqualMatrix(matrix3, matrix4, MatrixSize));
}

TEST(Strassen_OMP, speed_strassen_omp_16x16) {
  int MatrixSize = 16;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> matrix3(MatrixResultSize);
  std::vector<double> matrix4(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = i;
  }

  double time1, time2;
  time1 = omp_get_wtime();
  matrix3 = StrassenAlgorithm(matrix1, matrix2, MatrixSize);
  time2 = omp_get_wtime();
  time2 -= time1;
  std::cout << "Time of StrassenAlgorithm is " << time2 << std::endl;

  double time_omp1, time_omp2;
  time_omp1 = omp_get_wtime();
  StrassenAlgorithm_omp(matrix1, matrix2, &matrix4);
  time_omp2 = omp_get_wtime();
  time_omp2 -= time_omp1;
  std::cout << "Time of StrassenAlgorithm_omp is " << time_omp2 << std::endl;
  std::cout << "StrassenAlgorithm_omp is " << time2 / time_omp2 << " faster"
            << std::endl;

  ASSERT_TRUE(CheckEqualMatrix(matrix3, matrix4, MatrixSize));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
