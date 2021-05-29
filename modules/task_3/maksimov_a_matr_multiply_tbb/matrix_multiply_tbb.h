// Copyright 2021 Maksimov Andrey
#ifndef MODULES_TASK_3_MAKSIMOV_A_MATR_MULTIPLY_TBB_MATRIX_MULTIPLY_TBB_H_
#define MODULES_TASK_3_MAKSIMOV_A_MATR_MULTIPLY_TBB_MATRIX_MULTIPLY_TBB_H_
#include <vector>

std::vector<double> generateMatrix(int rows);
std::vector<double> mupltiplyMatrixByMatrix(
    std::vector<double> matr1,
    std::vector<double> matr2,
    int matr1Rows, int matr2Rows);
std::vector<double> mupltiplyMatrixByMatrixTbb(
    std::vector<double> matr1,
    std::vector<double> matr2,
    int matr1Rows, int matr2Rows,
    int threadsNum);
#endif  // MODULES_TASK_3_MAKSIMOV_A_MATR_MULTIPLY_TBB_MATRIX_MULTIPLY_TBB_H_
