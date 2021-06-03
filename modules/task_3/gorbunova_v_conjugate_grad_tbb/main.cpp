// Copyright 2021 Gorbunova Valeria
#include <stdio.h>
#include <gtest/gtest.h>
#include <tbb/blocked_range.h>
#include <tbb/tbb.h>
#include "./conjugate_gradient.h"



TEST(Conj_Grad, test_1_Rand_Matrix_10x10_4_proc) {
    int S = 10;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_gradP(A, B, S, 4));
}

TEST(Conj_Grad, test_2_Rand_Matrix_500x500_4_proc) {
    int S = 500;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_gradP(A, B, S, 4));
}

TEST(Conj_Grad, test_3_Rand_Matrix_500x500_10_proc) {
    int S = 500;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_gradP(A, B, S, 10));
}

TEST(Conj_Grad, test_4_Rand_Matrix_1000x1000_4_proc) {
    int S = 1000;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_gradP(A, B, S, 4));
}

TEST(Conj_Grad, test_5_Rand_Matrix_1000x1000_10_proc) {
    int S = 1000;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_gradP(A, B, S, 10));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
