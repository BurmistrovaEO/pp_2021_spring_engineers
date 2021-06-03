// Copyright 2021 Gorbunova Valeria
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include "../../../modules/task_2/gorbunova_v_conjugate_grad_omp/conjugate_gradient.h"

void Random_Matrix_A(int S, double** A) {
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < S; j++) {
            A[i][i] = std::rand() % 2;
            if ((i != j) && (i < j)) {
                A[i][j] = std::rand() % 100 - std::rand() % 100;
            }
            A[j][i] = A[i][j];
        }
    }
}

void Random_Matrix_B(int S, double* B) {
    for (int i = 0; i < S; i++) {
        B[i] = std::rand() % 100 - std::rand() % 100;
    }
}

void MatrixMultiplicate(double** A, int S, double* Ax, double* x1) {
#pragma omp parallel for
    for (int i = 0; i < S; i++) {
        Ax[i] = 0;
        for (int j = 0; j < S; j++)
            Ax[i] += A[i][j] * x1[j];
    }
}

double multiplicate(int S, double* x1, double* x2) {
    double result = 0;
#pragma omp parallel for reduction(+:result)
    for (int i = 0; i < S; i++) {
        result += x1[i] * x2[i];
    }
    return result;
}

double conj_grad(double** A, double* B, int S, int proc) {
    double* x = new double[S];
    double* d = new double[S];
    double* r = new double[S];
    double* Ax = new double[S];
    double alpha, beta, sumB, num, R, den;
    int i, IterNum = 1;

    double MaxIteration = 100000;
    for (i = 0; i < S; i++) {
        x[i] = 1;
    }
    omp_set_num_threads(proc);
    MatrixMultiplicate(A, S, Ax, x);
    for (int i = 0; i < S; i++) {
        r[i] = B[i] - Ax[i];
        d[i] = r[i];
    }
    for (sumB = 0, i = 0; i < S; i++) {
        sumB += B[i] * B[i];
    }

    int Iteration = 0;
    do {
        Iteration++;
        MatrixMultiplicate(A, S, Ax, d);
        num = multiplicate(S, r, r);
        den = multiplicate(S, Ax, d);
        alpha = num / den;
        for (i = 0; i < S; i++) {
            x[i] += alpha * d[i];
            r[i] -= alpha * Ax[i];
        }
        IterNum++;
        R = multiplicate(S, r, r);
        beta = R / num;
        for (i = 0; i < S; i++)
            d[i] = r[i] + beta * d[i];
    } while (R / sumB > 0.01 * 0.01 && Iteration < MaxIteration);

    int res = 0;
    double* A1 = new double[S];
    int F = S / 10;
    if (S / 10 == 0) {
        F = 1;
    }
    double fault = 0.3 * F;
    for (i = 0; i < S; i++) {
        A1[i] = 0;
        for (int j = 0; j < S; j++) {
            A1[i] += A[i][j]*x[j];
        }
        if (((abs(A1[i]) > abs(B[i])) && (abs(A1[i]) - abs(B[i]) > fault))||
        ((abs(A1[i]) <= abs(B[i])) && (abs(B[i]) - abs(A1[i]) > fault))) {
            res = 1;
            break;
        }
    }
    delete[] Ax;
    for (int i = 0; i < S; i++)
        delete[] A[i];
    delete[] A;
    delete[] A1;
    delete[] B;
    delete[] d;
    delete[] r;
    delete[] x;
    return res;
}
