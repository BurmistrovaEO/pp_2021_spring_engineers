// Copyright 2021 Gorbunova Valeria
#include <stdio.h>
#include <gtest/gtest.h>
#include <omp.h>
#include "./conjugate_gradient.h"

TEST(Conj_Grad, test_2_Rand_Matrix_10x10_10_proc) {
    int S = 10;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_grad(A, B, S, 10));
}

TEST(Conj_Grad, test_1_Rand_Matrix_5x5_5_proc) {
    int S = 5;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < S; j++) {
            std::cout << A[i][j] << "  ";
        }
        std::cout <<"  " << B[i] << std::endl;
    }
    ASSERT_EQ(0, conj_grad(A, B, S, 3));
}

TEST(Conj_Grad, test_3_Rand_Matrix_50x50_4_proc) {
    int S = 50;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_grad(A, B, S, 4));
}

TEST(Conj_Grad, test_4_Rand_Matrix_50x50_9_proc) {
    int S = 50;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_grad(A, B, S, 9));
}

TEST(Conj_Grad, test_5_Rand_Matrix_100x100_13_proc) {
    int S = 100;
    double** A = new double* [S];
    for (int count = 0; count < S; count++)
        A[count] = new double[S];
    Random_Matrix_A(S, A);
    double* B = new double[S];
    Random_Matrix_B(S, B);
    ASSERT_EQ(0, conj_grad(A, B, S, 13));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
