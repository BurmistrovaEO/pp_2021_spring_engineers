// Copyright 2021 Polina Kren
#ifndef MODULES_TASK_3_KREN_P_MATRIX_TBB_MATRIX_H_
#define MODULES_TASK_3_KREN_P_MATRIX_TBB_MATRIX_H_

#include <tbb/task_group.h>
#include <tbb/tbb.h>

#include <vector>

std::vector<double> AdditionMatrix(const std::vector<double>& matrix1,
                                   const std::vector<double>& matrix2,
                                   const int MatrixSize);

void AdditionMatrix(const std::vector<double>& matrix1,
                    const std::vector<double>& matrix2,
                    std::vector<double>* matrix3);

std::vector<double> SubtractionMatrix(const std::vector<double>& matrix1,
                                      const std::vector<double>& matrix2,
                                      const int MatrixSize);

void SubtractionMatrix(const std::vector<double>& matrix1,
                       const std::vector<double>& matrix2,
                       std::vector<double>* matrix3);

std::vector<double> MultiplicationMatrix(const std::vector<double>& matrix1,
                                         const std::vector<double>& matrix2,
                                         const int MatrixSize);

void MultiplicationMatrix(const std::vector<double>& matrix1,
                          const std::vector<double>& matrix2,
                          std::vector<double>* matrix3, const int MatrixSize);

bool CheckEqualMatrix(const std::vector<double>& matrix1,
                      const std::vector<double>& matrix2, const int MatrixSize);

int CheckMatrixSize(const int MatrixSize);

std::vector<double> ResizeMatrix(const std::vector<double>& matrix,
                                 const int MatrixSize);

void StrassenAlgorithm(const std::vector<double>& matrix1,
                       const std::vector<double>& matrix2,
                       std::vector<double>* result);
void SplitMatrixIntoPieces(const std::vector<double>& matrix,
                           std::vector<double>* A11, std::vector<double>* A22,
                           std::vector<double>* A12, std::vector<double>* A21,
                           const int MatrixSize);

void RebuildMatrix(std::vector<double>* MatrixResult,
                   const std::vector<double>& A11,
                   const std::vector<double>& A22,
                   const std::vector<double>& A12,
                   const std::vector<double>& A21, const int MatrixSize);

void StrassenAlgorithm_tbb(const std::vector<double>& matrix1,
                           const std::vector<double>& matrix2,
                           std::vector<double>* matrix3);

#endif  // MODULES_TASK_3_KREN_P_MATRIX_TBB_MATRIX_H_
