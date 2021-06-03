// Copyright 2021 Maksimov Andrey
#include <tbb/tbb.h>
#include <time.h>
#include <math.h>
#include <random>
#include <vector>
#include "../../modules/task_3/maksimov_a_matr_multiply_tbb/matrix_multiply_tbb.h"

std::vector<double> generateMatrix(int rows) {
    std::mt19937 gen;
    gen.seed(time(0));

    std::vector<double> matr(rows * rows);
    for (int i = 0; i < rows * rows; i++) {
        matr[i] = static_cast<double>(gen() / (RAND_MAX / 3)) / 2 + 0.5;
    }
    return matr;
}

std::vector<double> mupltiplyMatrixByMatrix(
    std::vector<double> matr1,
    std::vector<double> matr2,
    int matr1Rows, int matr2Rows) {
    if (matr1Rows != matr2Rows)
        throw "matr1Rows != matr2Rows";

    std::vector<double> result(matr1Rows * matr1Rows);
    for (int i = 0; i < matr1Rows; i++) {
        for (int j = 0; j < matr1Rows; j++) {
            result[i * matr1Rows + j] = 0;
            for (int k = 0; k < matr1Rows; k++) {
                result[i * matr1Rows + j]
                    += matr1[i * matr1Rows + k] * matr2[k * matr1Rows + j];
            }
        }
    }
    return result;
}

std::vector<double> mupltiplyMatrixByMatrixTbb(
    std::vector<double> matr1,
    std::vector<double> matr2,
    int matr1Rows, int matr2Rows,
    int threadsNum) {
    if (matr1Rows != matr2Rows)
        throw "matr1Rows != matr2Rows";
    int threadsNumSqrt = static_cast<int>(std::sqrt(threadsNum));
    if (threadsNumSqrt * threadsNumSqrt != threadsNum)
        throw "The number of threads is not a perfect square";

    int elems = matr1Rows / threadsNumSqrt;
    if (static_cast<int>(static_cast<double>(matr1Rows) / threadsNumSqrt) != elems)
        throw "Cannot be divided into blocks";

    std::vector<double> result(matr1Rows * matr1Rows);
    for (int i = 0; i < matr1Rows * matr1Rows; i++)
        result[i] = 0;

    tbb::task_scheduler_init init(threadsNum);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, threadsNum, 1),
        [&result, &threadsNumSqrt, &elems, &matr1Rows, &matr1, &matr2](const tbb::blocked_range<size_t>& range) {
        for (size_t cur = range.begin(); cur != range.end(); ++cur) {
            int iter;
            int i, j, k;
            int X, Y;
            int AblockX, AblockY, BblockX, BblockY;
            int threadId = cur;
            X = threadId % threadsNumSqrt;
            Y = threadId / threadsNumSqrt;

            AblockX = Y + X;
            if (AblockX >= threadsNumSqrt)
                AblockX -= threadsNumSqrt;
            AblockY = Y;

            BblockX = X;
            BblockY = Y + X;
            if (BblockY >= threadsNumSqrt)
                BblockY -= threadsNumSqrt;

            for (iter = 0; iter < threadsNumSqrt; iter++) {
                for (i = 0; i < elems; i++) {
                    for (j = 0; j < elems; j++) {
                        for (k = 0; k < elems; k++) {
                            double f = matr1[(i + AblockY * elems) * matr1Rows + k + AblockX * elems];
                            double s = matr2[(k + BblockY * elems) * matr1Rows + j + BblockX * elems];
                            result[(i + Y * elems) * matr1Rows + j + X * elems]
                                += f * s;
                        }
                    }
                }

                AblockX++;
                if (AblockX == threadsNumSqrt)
                    AblockX = 0;
                BblockY++;
                if (BblockY == threadsNumSqrt)
                    BblockY = 0;
            }
        }
    });
    init.terminate();
    return result;
}
