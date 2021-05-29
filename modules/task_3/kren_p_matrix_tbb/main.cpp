// Copyright 2021 Polina Kren
#include <time.h>

#include <vector>

#include "../../modules/task_3/kren_p_matrix_tbb/matrix.h"
#include "gtest/gtest.h"

TEST(Strassen_tbb, correct_strassen_algorithm_4x4) {
  int MatrixSize = 4;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);
  std::vector<double> MatrixResultTbb(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = 2;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = 1;
  }
  MatrixResult = MultiplicationMatrix(matrix1, matrix2, MatrixSize);
  StrassenAlgorithm_tbb(matrix1, matrix2, &MatrixResultTbb);
  ASSERT_TRUE(CheckEqualMatrix(MatrixResult, MatrixResultTbb, MatrixSize));
}

TEST(Strassen_tbb, correct_strassen_algorithm_8x8) {
  int MatrixSize = 8;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);
  std::vector<double> MatrixResultTbb(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = MatrixResultSize - i;
  }
  MatrixResult = MultiplicationMatrix(matrix1, matrix2, MatrixSize);
  StrassenAlgorithm_tbb(matrix1, matrix2, &MatrixResultTbb);
  ASSERT_TRUE(CheckEqualMatrix(MatrixResult, MatrixResultTbb, MatrixSize));
}

TEST(Strassen_tbb, correct_strassen_algorithm_16x16) {
  int MatrixSize = 16;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);
  std::vector<double> MatrixResultTbb(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = MatrixResultSize - i;
  }
  MatrixResult = MultiplicationMatrix(matrix1, matrix2, MatrixSize);
  StrassenAlgorithm_tbb(matrix1, matrix2, &MatrixResultTbb);
  ASSERT_TRUE(CheckEqualMatrix(MatrixResult, MatrixResultTbb, MatrixSize));
}

TEST(Strassen_tbb, correct_strassen_algorithm_16x16_1) {
  int MatrixSize = 16;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);
  std::vector<double> MatrixResultTbb(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = 2;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = 1;
  }
  MatrixResult = MultiplicationMatrix(matrix1, matrix2, MatrixSize);
  StrassenAlgorithm_tbb(matrix1, matrix2, &MatrixResultTbb);
  ASSERT_TRUE(CheckEqualMatrix(MatrixResult, MatrixResultTbb, MatrixSize));
}

TEST(Strassen_tbb, correct_strassen_algorithm_256x256) {
  int MatrixSize = 16;
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> matrix1(MatrixResultSize);
  std::vector<double> matrix2(MatrixResultSize);
  std::vector<double> MatrixResult(MatrixResultSize);
  std::vector<double> MatrixResultTbb(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix1[i] = i;
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    matrix2[i] = MatrixResultSize - i;
  }

  double time1, time2;
  time1 = clock();
  StrassenAlgorithm_tbb(matrix1, matrix2, &MatrixResultTbb);
  time2 = clock();
  time2 -= time1;
  time2 /= CLOCKS_PER_SEC;
  std::cout << " Time of StrassenAlgorithm is " << time2 << std::endl;

  double time_tbb1, time_tbb2;
  time_tbb1 = clock();
  StrassenAlgorithm(matrix1, matrix2, &MatrixResult);
  time_tbb2 = clock();
  time_tbb2 -= time_tbb1;
  time_tbb2 /= CLOCKS_PER_SEC;
  std::cout << " Time of StrassenAlgorithm_tbb is " << time_tbb2 << std::endl;
  std::cout << "StrassenAlgorithm_tbb in " << time2 / time_tbb2 << " faster"
            << std::endl;
  ASSERT_TRUE(CheckEqualMatrix(MatrixResult, MatrixResultTbb, MatrixSize));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
