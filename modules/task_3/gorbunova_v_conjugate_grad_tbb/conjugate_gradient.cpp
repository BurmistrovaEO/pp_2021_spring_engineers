/// Copyright 2021 Gorbunova Valeria
#include <stdio.h>
#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <stdlib.h>
#include <iostream>
#include "tbb/blocked_range2d.h"
#include "../../../modules/task_3/gorbunova_v_conjugate_grad_tbb/conjugate_gradient.h"


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

void MatrixMultiplicateP(double** A, int S, double* Ax, double* x1, int grainsize) {
    tbb::parallel_for(tbb::blocked_range2d<int>(0 , S, grainsize , 0, S, S),
        [&](tbb::blocked_range2d<int> r) {
            for (int i = r.rows().begin(); i < r.rows().end(); i++) {
                Ax[i] = 0;
                for (int j = r.cols().begin(); j < r.cols().end(); j++)
                    Ax[i] += A[i][j] * x1[j];
            }
        });
}

double multiplicate(int S, double* x1, double* x2) {
    double result = 0;
    for (int i = 0; i < S; i++) {
        result += x1[i] * x2[i];
    }
    return result;
}

double conj_gradP(double** A, double* B, int S, int proc) {
    double* x = new double[S];
    double* d = new double[S];
    double* r = new double[S];
    double* Ax = new double[S];
    double alpha, beta, sumB, num, R, den;
    int i, IterNum = 1;
    int grainsize = S / proc;
    tbb::task_scheduler_init init(proc);
    tbb::tick_count start2 = tbb::tick_count::now();

    double MaxIteration = 100000;
    for (i = 0; i < S; i++) {
        x[i] = std::rand() % 100 - std::rand() % 100;;
    }

    MatrixMultiplicateP(A, S, Ax, x, grainsize);
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
        MatrixMultiplicateP(A, S, Ax, d, grainsize);
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
    double fault = 0.2 * F;
    for (i = 0; i < S; i++) {
        A1[i] = 0;
        for (int j = 0; j < S; j++) {
            A1[i] += A[i][j] * x[j];
        }
        if (((abs(A1[i]) > abs(B[i])) && (abs(A1[i]) - abs(B[i]) > fault)) ||
            ((abs(A1[i]) <= abs(B[i])) && (abs(B[i]) - abs(A1[i]) > fault))) {
            res = 1;
            break;
        }
    }
    tbb::tick_count end2 = tbb::tick_count::now();
    std::cout << "Tbb time: " << (end2 - start2).seconds() << std::endl;
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
