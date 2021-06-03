// Copyright 2021 Kren Polina
#ifndef MODULES_TASK_2_KREN_P_MATRIX_OMP_MATRIX_H_
#define MODULES_TASK_2_KREN_P_MATRIX_OMP_MATRIX_H_

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

std::vector<double> StrassenAlgorithm(const std::vector<double>& matrix1,
                                      const std::vector<double>& matrix2,
                                      int MatrixSize);

void StrassenAlgorithm_omp(const std::vector<double>& matrix1,
                           const std::vector<double>& matrix2,
                           std::vector<double>* MatrixResult);

void SplitMatrixIntoPieces(const std::vector<double>& matrix,
                           std::vector<double>* A11, std::vector<double>* A22,
                           std::vector<double>* A12, std::vector<double>* A21,
                           const int MatrixSize);

void SplitMatrixIntoPiecesOmp(const std::vector<double>& matrix,
                              std::vector<double>* A11,
                              std::vector<double>* A22,
                              std::vector<double>* A12,
                              std::vector<double>* A21, const int MatrixSize);

std::vector<double> RebuildMatrix(const std::vector<double>& A11,
                                  const std::vector<double>& A22,
                                  const std::vector<double>& A12,
                                  const std::vector<double>& A21,
                                  const int SizeOfPiece);

void RebuildMatrixOmp(std::vector<double>* MatrixResult,
                      const std::vector<double>& A11,
                      const std::vector<double>& A22,
                      const std::vector<double>& A12,
                      const std::vector<double>& A21, const int MatrixSize);

#endif  // MODULES_TASK_2_KREN_P_MATRIX_OMP_MATRIX_H_
