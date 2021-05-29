// Copyright 2021 Zlobin George
#include <gtest/gtest.h>
#include <iostream>
#include "./matrix_ccs_mult.h"

// Tests for Constructor MatrixCCS(int nCollumns, int nRows, int nNotZero);
TEST(MatrixCCSConstructorFromNCollumnsNRowsNNotZero, ThrowIfIncorrectNCollumns) {
    ASSERT_ANY_THROW(MatrixCCS(-10, 10, 3));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsNNotZero, ThrowIfIncorrectNRows) {
    ASSERT_ANY_THROW(MatrixCCS(10, 0, 3));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsNNotZero, ThrowIfIncorrectNNotZero) {
    ASSERT_ANY_THROW(MatrixCCS(10, 10, -3));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsNNotZero, ThrowIfNNotZeroMoreThanAvailable) {
    ASSERT_ANY_THROW(MatrixCCS(10, 10, 101));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsNNotZero, NoThrowWithCorrectParameters) {
    ASSERT_NO_THROW(MatrixCCS(10, 10, 3));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsNNotZero, CanCreateMatrixWithZeroElements) {
    ASSERT_NO_THROW(MatrixCCS(10, 10, 0));
}

// Tests for Constructor MatrixCCS(nCollumns, nRows, matrix);
TEST(MatrixCCSConstructorFromNCollumnsNRowsMatrix, ThrowIfIncorrectNCollumns) {
    ASSERT_ANY_THROW(MatrixCCS(0, 2,
        { {2, 3}, {0, 0},
        {0, -1}, {-7, 0} }));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsMatrix, ThrowIfIncorrectNRows) {
    ASSERT_ANY_THROW(MatrixCCS(2, -3,
        { {0, 0}, {2, 3},
        {0, -1}, {-7, 0} }));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsMatrix, ThrowIfMatrixHasIncorrectSize) {
    ASSERT_ANY_THROW(MatrixCCS(2, 2,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} }));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsMatrix, NoThrowWithCorrectParameters) {
    ASSERT_NO_THROW(MatrixCCS(3, 3,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} }));
}

TEST(MatrixCCSConstructorFromNCollumnsNRowsMatrix, CanCreateMatrixWithZeroElements) {
    ASSERT_NO_THROW(MatrixCCS(3, 3,
        { {0, 0}, {0, 0}, {0, 0},
        {0, 0}, {0, 0}, {0, 0},
        {0, 0}, {0, 0}, {0, 0} }));
}

// Tests for Constructor MatrixCCS(int nCollumns, int nRows,
//                                 std::vector<int> collumnsIndexes,
//                                 std::vector<int> rows,
//                                 std::vector<std::complex<int>> values);
TEST(MatrixCCSConstructorFromDirectGivenFields, ThrowIfIncorrectNCollumns) {
    ASSERT_ANY_THROW(MatrixCCS(-2, 2,
        { 0, 0, 2 },
        { 0, 1 },
        { {-7, 0}, {0, -10} }));
}

TEST(MatrixCCSConstructorFromDirectGivenFields, ThrowIfIncorrectNRows) {
    ASSERT_ANY_THROW(MatrixCCS(2, 0,
        { 0, 0, 2 },
        { 0, 1 },
        { {-7, 0}, {0, -10} }));
}

TEST(MatrixCCSConstructorFromDirectGivenFields, ThrowIfGivenMoreValuesThanMatrixSize) {
    ASSERT_ANY_THROW(MatrixCCS(2, 0,
        { 0, 0, 2 },
        { 0, 1 },
        { {-7, 0}, {0, -10} }));
}

TEST(MatrixCCSConstructorFromDirectGivenFields, ThrowIfIncorrectCollumnsIndexes) {
    ASSERT_ANY_THROW(MatrixCCS(2, 2,
        { 0, 0, 2, 2 },
        { 0, 1 },
        { {-7, 0}, {0, -10} }));
}

TEST(MatrixCCSConstructorFromDirectGivenFields, ThrowIfRowsHasIncorrectValue) {
    ASSERT_ANY_THROW(MatrixCCS(2, 2,
        { 0, 0, 2 },
        { 0, 2 },
        { {-7, 0}, {0, -10} }));
}

TEST(MatrixCCSConstructorFromDirectGivenFields, NoThrowWithCorrectParameters) {
    ASSERT_NO_THROW(MatrixCCS(2, 2,
        { 0, 0, 2 },
        { 0, 1 },
        { {-7, 0}, {0, -10} }));
}

TEST(MatrixCCSConstructorFromDirectGivenFields, CanCreateMatrixWithZeroElements) {
    ASSERT_NO_THROW(MatrixCCS(2, 2,
        { 0, 0, 0},
        {},
        {}));
}

// Tests for constructors MatrixCCS are equivalent
TEST(MatrixCCSConstructorAreEquivalent, MatrixWithZeroElements) {
    EXPECT_EQ(MatrixCCS(3, 3,
        { {0, 0}, {0, 0}, {0, 0},
        {0, 0}, {0, 0}, {0, 0},
        {0, 0}, {0, 0}, {0, 0} }),
              MatrixCCS(3, 3,
        { 0, 0, 0, 0},
        {},
        {}));
}

TEST(MatrixCCSConstructorAreEquivalent, SquareMatrix) {
    EXPECT_EQ(MatrixCCS(3, 3,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} }),
              MatrixCCS(3, 3,
        { 0, 1, 3, 5 },
        { 1, 0, 1, 0, 2 },
        { {0, -1}, {2, 3}, {-7, 0},
        {7, 0}, {-3, 10} }));
}

TEST(MatrixCCSConstructorAreEquivalent, Row) {
    EXPECT_EQ(MatrixCCS(3, 1,
        { {0, 0}, {2, 3}, {7, 0} }),
              MatrixCCS(3, 1,
        { 0, 0, 1, 2 },
        { 0, 0 },
        { {2, 3}, {7, 0} }));
}

TEST(MatrixCCSConstructorAreEquivalent, Collumn) {
    EXPECT_EQ(MatrixCCS(1, 3,
        { {0, 0},
        {0, -1},
        {0, 0} }),
              MatrixCCS(1, 3,
        { 0, 1 },
        { 1 },
        { {0, -1} }));
}

TEST(MatrixCCSConstructorAreEquivalent, NotSquareMatrix) {
    EXPECT_EQ(MatrixCCS(3, 2,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0} }),
              MatrixCCS(3, 2,
        { 0, 1, 3, 4 },
        { 1, 0, 1, 0 },
        { {0, -1}, {2, 3}, {-7, 0},
        {7, 0} }));
}

TEST(MatrixCCSConstructorAreEquivalent, DifferentSizesNotSame) {
    EXPECT_FALSE(MatrixCCS(3, 3,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} }) ==
              MatrixCCS(3, 2,
        { 0, 1, 3, 4 },
        { 1, 0, 1, 0 },
        { {0, -1}, {2, 3}, {-7, 0},
        {7, 0} }));
}

TEST(MatrixCCSConstructorAreEquivalent, DifferentValuesNotSame) {
    EXPECT_FALSE(MatrixCCS(3, 3,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} }) ==
              MatrixCCS(3, 3,
        { 0, 1, 3, 5 },
        { 1, 0, 1, 0, 2 },
        { {0, -1}, {2, 3}, {-7, 0},
        {7, 0}, {-3, 10} }));
}

// Tests for FillRandom(unsigned seed = 132, int min = -1000, int max = 1000);
TEST(MatrixCCSFillWithRandomValues, CanFillMatrix) {
    MatrixCCS A = MatrixCCS(3, 3, 0);
    ASSERT_NO_THROW(A.FillRandom());
}

TEST(MatrixCCSFillWithRandomValues, SameSeedSameMatrix) {
    MatrixCCS A = MatrixCCS(3, 3, 0);
    MatrixCCS B = A;
    A.FillRandom();
    B.FillRandom();
    EXPECT_EQ(A, B);
}

// Tests for Transpose()
TEST(MatrixCCSTranspose, CanTranspose1) {
    MatrixCCS A = MatrixCCS(3, 3,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} });
    ASSERT_NO_THROW(A.Transpose());
}
TEST(MatrixCCSTranspose, CanTranspose2) {
    MatrixCCS A = MatrixCCS(10, 5, 30);
    A.FillRandom();
    ASSERT_NO_THROW(A.Transpose());
}
TEST(MatrixCCSTranspose, SquareMatrix) {
    MatrixCCS A = MatrixCCS(3, 3,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {-3, 10} });
    MatrixCCS AT = A.Transpose();
    MatrixCCS ATref = MatrixCCS(3, 3,
        { {0, 0}, {0, -1}, {0, 0},
        {2, 3}, {-7, 0}, {0, 0},
        {7, 0}, {0, 0}, {-3, 10} });
    EXPECT_EQ(ATref, AT);
    EXPECT_EQ(A, AT.Transpose());
}

TEST(MatrixCCSTranspose, Row) {
    MatrixCCS A = MatrixCCS(3, 1,
        { {0, 0}, {2, 3}, {7, 0} });
    MatrixCCS AT = A.Transpose();
    MatrixCCS ATref = MatrixCCS(1, 3,
        { {0, 0},
        {2, 3},
        {7, 0} });
    EXPECT_EQ(ATref, AT);
    EXPECT_EQ(A, AT.Transpose());
}

TEST(MatrixCCSTranspose, NotSquareMatrix) {
    MatrixCCS A = MatrixCCS(3, 2,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0} });
    MatrixCCS AT = A.Transpose();
    MatrixCCS ATref = MatrixCCS(2, 3,
        { {0, 0}, {0, -1},
        {2, 3}, {-7, 0},
        {7, 0}, {0, 0} });
    EXPECT_EQ(ATref, AT);
    EXPECT_EQ(A, AT.Transpose());
}

// Tests for Multiplication
TEST(MatrixCCSMultiplication, ThrowIfNotCorrespondingMatrices) {
    MatrixCCS A = MatrixCCS(3, 2,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0} });
    MatrixCCS B = MatrixCCS(2, 4,
        { {0, 0}, {0, -1},
        {2, 3}, {-7, 0},
        {7, 0}, {0, 0},
        {0, -10}, {20, 0} });
    ASSERT_ANY_THROW(TBBMultiplication(A, B));
}

TEST(MatrixCCSMultiplication, CanMultiplyCorrespondingMatrices) {
    MatrixCCS A = MatrixCCS(3, 2,
        { {0, 0}, {2, 3}, {7, 0},
        {0, -1}, {-7, 0}, {0, 0} });
    MatrixCCS B = MatrixCCS(2, 3,
        { {0, 0}, {0, -1},
        {2, 3}, {-7, 0},
        {7, 0}, {0, 0} });
    ASSERT_NO_THROW(TBBMultiplication(A, B));
}

TEST(MatrixCCSMultiplication, MyltiplyNumbers) {
    MatrixCCS A = MatrixCCS(1, 1,
        { {20, -30} });
    MatrixCCS B = MatrixCCS(1, 1,
        { {86, 360} });
    MatrixCCS AB = MatrixCCS(1, 1,
        { {12520, 4620} });
    ASSERT_EQ(AB, TBBMultiplication(A, B));
}

TEST(MatrixCCSMultiplication, MyltiplyOnZeroMatrix) {
    MatrixCCS A = MatrixCCS(3, 2,
        { {0, 0}, {20, 0}, {7, 0},
        {8, 0}, {-7, 0}, {0, 0} });
    MatrixCCS B = MatrixCCS(2, 3,
        { {0, 0}, {0, 0},
        {0, 0}, {0, 0},
        {0, 0}, {0, 0} });
    MatrixCCS AB = MatrixCCS(2, 2,
        { {0, 0}, {0, 0},
        {0, 0}, {0, 0} });
    ASSERT_EQ(AB, TBBMultiplication(A, B));
}

TEST(MatrixCCSMultiplication, MyltiplyOnSingularMatrix) {
    MatrixCCS A = MatrixCCS(2, 2,
        { {1, 0}, {0, 0},
        {0, 0}, {1, 0} });
    MatrixCCS B = MatrixCCS(2, 2,
        { {3, 20}, {4, 0},
        {0, -3}, {0, 0} });
    ASSERT_EQ(B, TBBMultiplication(A, B));
}

TEST(MatrixCCSMultiplication, MatrixMultiplicationExample1) {
    MatrixCCS A = MatrixCCS(3, 2,
        { {0, 0}, {20, 0}, {7, 0},
        {8, 0}, {-7, 0}, {0, 0} });
    MatrixCCS B = MatrixCCS(2, 3,
        { {0, 0}, {0, 0},
        {4, 0}, {-7, 0},
        {0, 0}, {0, 0} });
    MatrixCCS AB = MatrixCCS(2, 2,
        { {80, 0}, {-140, 0},
        {-28, 0}, {49, 0} });
    ASSERT_EQ(AB, TBBMultiplication(A, B));
}

TEST(MatrixCCSMultiplication, MatrixMultiplicationExample2) {
    MatrixCCS A = MatrixCCS(3, 4,
        { {0, 0}, {0, 0}, {0, 0},
        {0, 0}, {-7, 0}, {0, 0},
        {0, 0}, {0, 0}, {4, 7},
        {4, -3}, {0, 0}, {0, 0} });
    MatrixCCS B = MatrixCCS(2, 3,
        { {0, -3}, {0, 0},
        {0, 0}, {0, 0},
        {0, 0}, {7, -10} });
    MatrixCCS AB = MatrixCCS(2, 4,
        { {0, 0}, {0, 0},
        {0, 0}, {0, 0},
        {0, 0}, {98, 9},
        {-9, -12}, {0, 0} });
    ASSERT_EQ(AB, TBBMultiplication(A, B));
}

// Tests for compare with reference sequential realization
TEST(MatrixCCSCompareRealizations, CompareLeftRow) {
    MatrixCCS A = MatrixCCS(10, 1, 8);
    MatrixCCS B = MatrixCCS(5, 10, 30);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}
TEST(MatrixCCSCompareRealizations, CompareRightCollumn) {
    MatrixCCS A = MatrixCCS(10, 5, 30);
    MatrixCCS B = MatrixCCS(1, 10, 8);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}
TEST(MatrixCCSCompareRealizations, CompareLeftBigger) {
    MatrixCCS A = MatrixCCS(10, 20, 50);
    MatrixCCS B = MatrixCCS(5, 10, 30);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}
TEST(MatrixCCSCompareRealizations, CompareRightBigger) {
    MatrixCCS A = MatrixCCS(10, 5, 30);
    MatrixCCS B = MatrixCCS(20, 10, 50);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}
TEST(MatrixCCSCompareRealizations, CompareAlmostFull) {
    MatrixCCS A = MatrixCCS(10, 15, 149);
    MatrixCCS B = MatrixCCS(20, 10, 197);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}
TEST(MatrixCCSCompareRealizations, CompareFull) {
    MatrixCCS A = MatrixCCS(10, 20, 200);
    MatrixCCS B = MatrixCCS(15, 10, 150);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}
TEST(MatrixCCSCompareRealizations, CompareSparse) {
    MatrixCCS A = MatrixCCS(10, 20, 7);
    MatrixCCS B = MatrixCCS(15, 10, 5);
    A.FillRandom();
    B.FillRandom();
    ASSERT_EQ(SeqMultiplication(A, B), TBBMultiplication(A, B));
}

// Test for performance measure
// TEST(MatrixCCSPerformance, Sequential) {
//     MatrixCCS A = MatrixCCS(4000, 4000, 1000000);
//     MatrixCCS B = MatrixCCS(4000, 4000, 1000000);
//     A.FillRandom();
//     B.FillRandom();
//     tbb::tick_count start = tbb::tick_count::now();
//     ASSERT_NO_THROW(SeqMultiplication(A, B));
//     tbb::tick_count end = tbb::tick_count::now();
//     std::cout << "[     TIME ] " << (end - start).seconds() << "\tseconds\n";
// }
// TEST(MatrixCCSPerformance, TBB1Thread) {
//     MatrixCCS A = MatrixCCS(4000, 4000, 1000000);
//     MatrixCCS B = MatrixCCS(4000, 4000, 1000000);
//     A.FillRandom();
//     B.FillRandom();
//     tbb::task_scheduler_init init(1);
//     tbb::tick_count start = tbb::tick_count::now();
//     ASSERT_NO_THROW(TBBMultiplication(A, B));
//     tbb::tick_count end = tbb::tick_count::now();
//     init.terminate();
//     std::cout << "[     TIME ] " << (end - start).seconds() << "\tseconds\n";
// }
// TEST(MatrixCCSPerformance, TBB2Thread) {
//     MatrixCCS A = MatrixCCS(4000, 4000, 1000000);
//     MatrixCCS B = MatrixCCS(4000, 4000, 1000000);
//     A.FillRandom();
//     B.FillRandom();
//     tbb::task_scheduler_init init(2);
//     tbb::tick_count start = tbb::tick_count::now();
//     ASSERT_NO_THROW(TBBMultiplication(A, B));
//     tbb::tick_count end = tbb::tick_count::now();
//     init.terminate();
//     std::cout << "[     TIME ] " << (end - start).seconds() << "\tseconds\n";
// }
// TEST(MatrixCCSPerformance, TBB4Thread) {
//     MatrixCCS A = MatrixCCS(4000, 4000, 1000000);
//     MatrixCCS B = MatrixCCS(4000, 4000, 1000000);
//     A.FillRandom();
//     B.FillRandom();
//     tbb::task_scheduler_init init(4);
//     tbb::tick_count start = tbb::tick_count::now();
//     ASSERT_NO_THROW(TBBMultiplication(A, B));
//     tbb::tick_count end = tbb::tick_count::now();
//     init.terminate();
//     std::cout << "[     TIME ] " << (end - start).seconds() << "\tseconds\n";
// }
// TEST(MatrixCCSPerformance, TBB9Thread) {
//     MatrixCCS A = MatrixCCS(4000, 4000, 1000000);
//     MatrixCCS B = MatrixCCS(4000, 4000, 1000000);
//     A.FillRandom();
//     B.FillRandom();
//     tbb::task_scheduler_init init(9);
//     tbb::tick_count start = tbb::tick_count::now();
//     ASSERT_NO_THROW(TBBMultiplication(A, B));
//     tbb::tick_count end = tbb::tick_count::now();
//     init.terminate();
//     std::cout << "[     TIME ] " << (end - start).seconds() << "\tseconds\n";
// }

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
