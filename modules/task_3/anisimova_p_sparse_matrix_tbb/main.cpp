// Copyright 2021 Anisimova Polina
#include <gtest/gtest.h>
#include <tbb/tbb.h>

#include <cmath>
#include <ctime>
#include <vector>

#include "../../../modules/task_3/anisimova_p_sparse_matrix_tbb/sparsematrix.h"

TEST(SparseMatrixTest, MatrixIsEqualToItself1) {
  SparseMatrix A = SparseMatrix({1, 1, 1}, {0, 1, 2, 3}, {0, 1, 2}, 3, 3);
  ASSERT_TRUE(A == A);
}

TEST(SparseMatrixTest, MatrixIsEqualToItself2) {
  SparseMatrix A =
      SparseMatrix({5, 8, 1, 2, 6}, {0, 2, 3, 5, 5}, {0, 3, 2, 1, 3}, 5, 4);
  ASSERT_TRUE(A == A);
}

TEST(SparseMatrixTest, NoThrowTranspose) {
  SparseMatrix A =
      SparseMatrix({1, 7, 1, 2, 6}, {0, 2, 3, 5, 5}, {0, 3, 2, 1, 3}, 5, 4);
  ASSERT_NO_THROW(A.Transpose());
}

TEST(SparseMatrixTest, MatrixMult1) {
  SparseMatrix A = SparseMatrix({1, 1, 1}, {0, 1, 2, 3}, {0, 1, 2}, 3, 3);
  SparseMatrix B = SparseMatrix({7, 8, 9}, {0, 1, 2, 3}, {0, 1, 2}, 3, 3);
  ASSERT_TRUE(B == A * B);
}

TEST(SparseMatrixTest, MatrixMult2) {
  SparseMatrix A =
      SparseMatrix({5, 8, 1, 5, 6}, {0, 2, 3, 5, 5}, {0, 3, 2, 1, 3}, 5, 4);
  SparseMatrix B =
      SparseMatrix({1, 1, 1, 1}, {0, 1, 2, 3, 4}, {0, 1, 2, 3}, 4, 4);
  ASSERT_TRUE(A == A * B);
}

TEST(SparseMatrixTest, MatrixMult3) {
  SparseMatrix A = SparseMatrix({1, 1}, {0, 1, 1, 2}, {0, 2}, 2, 3);
  SparseMatrix B = SparseMatrix({7, 8, 9}, {0, 1, 2, 3}, {0, 1, 2}, 3, 3);
  SparseMatrix C = SparseMatrix({7, 9}, {0, 1, 1, 2}, {0, 2}, 2, 3);
  SparseMatrix D = A * B;
  ASSERT_TRUE(C == D);
}

TEST(SparseMatrixTest, MatrixMultEqualsZero) {
  SparseMatrix B =
      SparseMatrix({1, 7, 1, 2, 6}, {0, 2, 3, 5, 5}, {2, 3, 2, 1, 3}, 5, 4);
  SparseMatrix A = SparseMatrix({1}, {0, 1, 1, 1, 1}, {0}, 1, 4);
  SparseMatrix C = SparseMatrix({}, {0, 0, 0, 0, 0}, {}, 0, 4);
  SparseMatrix D = A * B;
  ASSERT_TRUE(C == D);
}

TEST(SparseMatrixTest, MatrixRandomEqualsWithSameSeed1) {
  SparseMatrix A = SparseMatrix(17, 8, 35);
  SparseMatrix B = SparseMatrix(17, 8, 35);
  ASSERT_TRUE(A == B);
}

TEST(SparseMatrixTest, MatrixRandomMult) {
  SparseMatrix A = SparseMatrix(1000, 50);
  SparseMatrix B = SparseMatrix(1000, 50);
  ASSERT_NO_THROW(A * B);
}

TEST(SparseMatrixTest, MatrixRandomMultParallel1) {
  SparseMatrix A = SparseMatrix(11000, 5000);
  SparseMatrix B = SparseMatrix(11000, 5000);

  auto t1 = tbb::tick_count::now();
  SparseMatrix C1 = A * B;
  std::cout << (tbb::tick_count::now() - t1).seconds() << std::endl;

  t1 = tbb::tick_count::now();
  SparseMatrix C2 = A.ParallelMult(B, 4);
  std::cout << (tbb::tick_count::now() - t1).seconds() << std::endl;

  bool flag = true;
  for (size_t i = 0; i < C1.values.size(); i++) {
    if (fabs(C1.values[i] - C2.values[i]) >=
        std::numeric_limits<double>::epsilon()) {
      flag = false;
      break;
    }
  }

  ASSERT_TRUE(flag);
}

TEST(SparseMatrixTest, MatrixRandomMultParallel2) {
  SparseMatrix A = SparseMatrix(3000, 1001);
  SparseMatrix B = SparseMatrix(3000, 1001);

  auto t1 = tbb::tick_count::now();
  SparseMatrix C1 = A * B;
  std::cout << (tbb::tick_count::now() - t1).seconds() << std::endl;

  t1 = tbb::tick_count::now();
  SparseMatrix C2 = A.ParallelMult(B, 4);
  std::cout << (tbb::tick_count::now() - t1).seconds() << std::endl;

  bool flag = true;
  for (size_t i = 0; i < C1.values.size(); i++) {
    if (fabs(C1.values[i] - C2.values[i]) >=
        std::numeric_limits<double>::epsilon()) {
      flag = false;
      break;
    }
  }

  ASSERT_TRUE(flag);
}

TEST(SparseMatrixTest, MatrixRandomMultParallel3) {
  SparseMatrix A = SparseMatrix(10000, 4000);
  SparseMatrix B = SparseMatrix(10000, 4000);

  auto t1 = tbb::tick_count::now();
  SparseMatrix C1 = A * B;
  std::cout << (tbb::tick_count::now() - t1).seconds() << std::endl;

  t1 = tbb::tick_count::now();
  SparseMatrix C2 = A.ParallelMult(B, 4);
  std::cout << (tbb::tick_count::now() - t1).seconds() << std::endl;

  bool flag = true;
  for (size_t i = 0; i < C1.values.size(); i++) {
    if (fabs(C1.values[i] - C2.values[i]) >=
        std::numeric_limits<double>::epsilon()) {
      flag = false;
      break;
    }
  }

  ASSERT_TRUE(flag);
}

TEST(SparseMatrixTest, EmptyMatrixRandomMultParallel4) {
  SparseMatrix A = SparseMatrix(15000, 0);
  SparseMatrix B = SparseMatrix(15000, 0);

  SparseMatrix C1 = A * B;
  SparseMatrix C2 = A.ParallelMult(B, 4);

  ASSERT_TRUE(true == C2.values.empty());
}

TEST(SparseMatrixTest, EmptyMatrixRandomMultParallel5) {
  SparseMatrix A = SparseMatrix(10000, 0);
  SparseMatrix B = SparseMatrix(10000, 0);

  SparseMatrix C1 = A * B;
  SparseMatrix C2 = A.ParallelMult(B, 4);

  ASSERT_TRUE(true == C2.values.empty());
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
