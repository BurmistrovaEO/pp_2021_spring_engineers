// Copyright 2021 Kren Polina

#include "../../modules/task_2/kren_p_matrix_omp/matrix.h"

#include <math.h>
#include <omp.h>

#include <algorithm>
#include <vector>

std::vector<double> AdditionMatrix(const std::vector<double>& matrix1,
                                   const std::vector<double>& matrix2,
                                   const int MatrixSize) {
  int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> MatrixResult(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    MatrixResult[i] += matrix1[i];
  }
  for (int i = 0; i < MatrixResultSize; i++) {
    MatrixResult[i] += matrix2[i];
  }
  return MatrixResult;
}

void AdditionMatrix(const std::vector<double>& matrix1,
                    const std::vector<double>& matrix2,
                    std::vector<double>* matrix3) {
  for (unsigned int i = 0; i < matrix3->size(); i++) {
    matrix3->at(i) += matrix1[i];
  }
  for (unsigned int i = 0; i < matrix3->size(); i++) {
    matrix3->at(i) += matrix2[i];
  }
}

std::vector<double> SubtractionMatrix(const std::vector<double>& matrix1,
                                      const std::vector<double>& matrix2,
                                      const int MatrixSize) {
  const int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> MatrixResult(MatrixResultSize);
  for (int i = 0; i < MatrixResultSize; i++) {
    MatrixResult[i] = matrix1[i] - matrix2[i];
  }
  return MatrixResult;
}

void SubtractionMatrix(const std::vector<double>& matrix1,
                       const std::vector<double>& matrix2,
                       std::vector<double>* matrix3) {
  for (unsigned int i = 0; i < matrix3->size(); i++) {
    matrix3->at(i) += matrix1[i];
  }
  for (unsigned int i = 0; i < matrix3->size(); i++) {
    matrix3->at(i) -= matrix2[i];
  }
}

std::vector<double> MultiplicationMatrix(const std::vector<double>& matrix1,
                                         const std::vector<double>& matrix2,
                                         const int MatrixSize) {
  const int MatrixResultSize = MatrixSize * MatrixSize;
  std::vector<double> MatrixResult(MatrixResultSize);
  for (int i = 0; i < MatrixSize; i++) {
    for (int j = 0; j < MatrixSize; j++) {
      MatrixResult[i * MatrixSize + j] = 0;
      for (int k = 0; k < MatrixSize; k++) {
        MatrixResult[i * MatrixSize + j] +=
            matrix1[i * MatrixSize + k] * matrix2[k * MatrixSize + j];
      }
    }
  }
  return MatrixResult;
}

void MultiplicationMatrix(const std::vector<double>& matrix1,
                          const std::vector<double>& matrix2,
                          std::vector<double>* matrix3, const int MatrixSize) {
  for (int i = 0; i < MatrixSize; i++) {
    for (int j = 0; j < MatrixSize; j++) {
      for (int k = 0; k < MatrixSize; k++) {
        matrix3->at(i * MatrixSize + j) +=
            matrix1[i * MatrixSize + k] * matrix2[k * MatrixSize + j];
      }
    }
  }
}

bool CheckEqualMatrix(const std::vector<double>& matrix1,
                      const std::vector<double>& matrix2,
                      const int MatrixSize) {
  bool MatrixIsEqual = true;
  for (int i = 0; i < MatrixSize; i++)
    if (matrix1[i] != matrix2[i]) {
      MatrixIsEqual = false;
      break;
    }
  return MatrixIsEqual;
}

int CheckMatrixSize(const int MatrixSize) {
  int MatrixResultSize = 2;
  while (MatrixSize > MatrixResultSize) {
    MatrixResultSize = MatrixResultSize * 2;
  }
  return MatrixResultSize;
}

std::vector<double> ResizeMatrix(const std::vector<double>& matrix,
                                 const int MatrixSize) {
  const int MatrixResultSize = CheckMatrixSize(MatrixSize);
  std::vector<double> MatrixResult(MatrixResultSize * MatrixResultSize);
  for (int i = 0; i < MatrixResultSize * MatrixResultSize; i++) {
    MatrixResult[i] = 0;
  }
  for (int i = 0; i < MatrixSize; i++) {
    for (int j = 0; j < MatrixSize; j++) {
      MatrixResult[i * MatrixResultSize + j] = matrix[i * MatrixSize + j];
    }
  }
  return MatrixResult;
}

void SplitMatrixIntoPieces(const std::vector<double>& matrix,
                           std::vector<double>* A11, std::vector<double>* A22,
                           std::vector<double>* A12, std::vector<double>* A21,
                           const int MatrixSize) {
  const int SizeOfPiece = MatrixSize / 2;
  for (int i = 0; i < SizeOfPiece; i++) {
    for (int j = 0; j < SizeOfPiece; j++) {
      A11->at(i * SizeOfPiece + j) = matrix[i * MatrixSize + j];
      A12->at(i * SizeOfPiece + j) = matrix[i * MatrixSize + j + SizeOfPiece];
      A21->at(i * SizeOfPiece + j) = matrix[(i + SizeOfPiece) * MatrixSize + j];
      A22->at(i * SizeOfPiece + j) =
          matrix[(i + SizeOfPiece) * MatrixSize + j + SizeOfPiece];
    }
  }
}

void SplitMatrixIntoPiecesOmp(const std::vector<double>& matrix,
                              std::vector<double>* A11,
                              std::vector<double>* A22,
                              std::vector<double>* A12,
                              std::vector<double>* A21, const int MatrixSize) {
  const int SizeOfPiece = MatrixSize / 2;
#pragma omp parallel for default(shared)
  for (int i = 0; i < SizeOfPiece; i++) {
    for (int j = 0; j < SizeOfPiece; j++) {
      A11->at(i * SizeOfPiece + j) = matrix[i * MatrixSize + j];
      A12->at(i * SizeOfPiece + j) = matrix[i * MatrixSize + j + SizeOfPiece];
      A21->at(i * SizeOfPiece + j) = matrix[(i + SizeOfPiece) * MatrixSize + j];
      A22->at(i * SizeOfPiece + j) =
          matrix[(i + SizeOfPiece) * MatrixSize + j + SizeOfPiece];
    }
  }
}

std::vector<double> RebuildMatrix(const std::vector<double>& A11,
                                  const std::vector<double>& A22,
                                  const std::vector<double>& A12,
                                  const std::vector<double>& A21,
                                  const int MatrixSize) {
  const int MatrixResultSize = MatrixSize * 2;
  std::vector<double> MatrixResult(MatrixResultSize * MatrixResultSize);
  for (int i = 0; i < MatrixSize; i++) {
    for (int j = 0; j < MatrixSize; j++) {
      MatrixResult[i * MatrixResultSize + j] = A11[i * MatrixSize + j];
      MatrixResult[i * MatrixResultSize + j + MatrixSize] =
          A12[i * MatrixSize + j];
      MatrixResult[(i + MatrixSize) * MatrixResultSize + j] =
          A21[i * MatrixSize + j];
      MatrixResult[(i + MatrixSize) * MatrixResultSize + j + MatrixSize] =
          A22[i * MatrixSize + j];
    }
  }
  return MatrixResult;
}

void RebuildMatrixOmp(std::vector<double>* MatrixResult,
                      const std::vector<double>& A11,
                      const std::vector<double>& A22,
                      const std::vector<double>& A12,
                      const std::vector<double>& A21, const int MatrixSize) {
  const int MatrixResultSize = MatrixSize * 2;
#pragma omp parallel for default(shared)
  for (int i = 0; i < MatrixSize; i++) {
    for (int j = 0; j < MatrixSize; j++) {
      MatrixResult->at(i * MatrixResultSize + j) = A11[i * MatrixSize + j];
      MatrixResult->at(i * MatrixResultSize + j + MatrixSize) =
          A12[i * MatrixSize + j];
      MatrixResult->at((i + MatrixSize) * MatrixResultSize + j) =
          A21[i * MatrixSize + j];
      MatrixResult->at((i + MatrixSize) * MatrixResultSize + j + MatrixSize) =
          A22[i * MatrixSize + j];
    }
  }
}

std::vector<double> StrassenAlgorithm(const std::vector<double>& matrix1,
                                      const std::vector<double>& matrix2,
                                      int MatrixSize) {
  if (MatrixSize <= 2) {
    return MultiplicationMatrix(matrix1, matrix2, MatrixSize);
  }
  MatrixSize = MatrixSize / 2;
  int MatrixResultSize = MatrixSize * MatrixSize;

  std::vector<double> A11(MatrixResultSize);
  std::vector<double> A12(MatrixResultSize);
  std::vector<double> A21(MatrixResultSize);
  std::vector<double> A22(MatrixResultSize);
  std::vector<double> B11(MatrixResultSize);
  std::vector<double> B12(MatrixResultSize);
  std::vector<double> B21(MatrixResultSize);
  std::vector<double> B22(MatrixResultSize);

  SplitMatrixIntoPieces(matrix1, &A11, &A22, &A12, &A21, MatrixSize * 2);
  SplitMatrixIntoPieces(matrix2, &B11, &B22, &B12, &B21, MatrixSize * 2);

  std::vector<double> p1 =
      StrassenAlgorithm(AdditionMatrix(A11, A22, MatrixSize),
                        AdditionMatrix(B11, B22, MatrixSize), MatrixSize);
  std::vector<double> p2 =
      StrassenAlgorithm(AdditionMatrix(A21, A22, MatrixSize), B11, MatrixSize);
  std::vector<double> p3 = StrassenAlgorithm(
      A11, SubtractionMatrix(B12, B22, MatrixSize), MatrixSize);
  std::vector<double> p4 = StrassenAlgorithm(
      A22, SubtractionMatrix(B21, B11, MatrixSize), MatrixSize);
  std::vector<double> p5 =
      StrassenAlgorithm(AdditionMatrix(A11, A12, MatrixSize), B22, MatrixSize);
  std::vector<double> p6 =
      StrassenAlgorithm(SubtractionMatrix(A21, A11, MatrixSize),
                        AdditionMatrix(B11, B12, MatrixSize), MatrixSize);
  std::vector<double> p7 =
      StrassenAlgorithm(SubtractionMatrix(A12, A22, MatrixSize),
                        AdditionMatrix(B21, B22, MatrixSize), MatrixSize);

  std::vector<double> c11(MatrixResultSize);
  std::vector<double> c12(MatrixResultSize);
  std::vector<double> c21(MatrixResultSize);
  std::vector<double> c22(MatrixResultSize);

  c11 = AdditionMatrix(AdditionMatrix(p1, p4, MatrixSize),
                       SubtractionMatrix(p7, p5, MatrixSize), MatrixSize);
  c12 = AdditionMatrix(p3, p5, MatrixSize);
  c21 = AdditionMatrix(p2, p4, MatrixSize);
  c22 = AdditionMatrix(SubtractionMatrix(p1, p2, MatrixSize),
                       AdditionMatrix(p3, p6, MatrixSize), MatrixSize);

  return RebuildMatrix(c11, c22, c12, c21, MatrixSize);
}

void StrassenAlgorithm_omp(const std::vector<double>& matrix1,
                           const std::vector<double>& matrix2,
                           std::vector<double>* MatrixResult) {
  int MatrixSize = sqrt(MatrixResult->size());

  if (MatrixSize <= 2) {
    return MultiplicationMatrix(matrix1, matrix2, MatrixResult, MatrixSize);
  }

  MatrixSize = MatrixSize / 2;
  int MatrixResultSize = MatrixSize * MatrixSize;

  std::vector<double> a11(MatrixResultSize);
  std::vector<double> a21(MatrixResultSize);
  std::vector<double> a12(MatrixResultSize);
  std::vector<double> a22(MatrixResultSize);
  std::vector<double> b11(MatrixResultSize);
  std::vector<double> b12(MatrixResultSize);
  std::vector<double> b21(MatrixResultSize);
  std::vector<double> b22(MatrixResultSize);

  SplitMatrixIntoPiecesOmp(matrix1, &a11, &a22, &a12, &a21, MatrixSize * 2);
  SplitMatrixIntoPiecesOmp(matrix2, &b11, &b22, &b12, &b21, MatrixSize * 2);

  std::vector<double> p1(MatrixResultSize);
  std::vector<double> p2(MatrixResultSize);
  std::vector<double> p3(MatrixResultSize);
  std::vector<double> p4(MatrixResultSize);
  std::vector<double> p5(MatrixResultSize);
  std::vector<double> p6(MatrixResultSize);
  std::vector<double> p7(MatrixResultSize);

#pragma omp parallel sections default(shared)
  {
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      std::vector<double> MatrixTemp2(MatrixResultSize);
      AdditionMatrix(a11, a22, &MatrixTemp1);
      AdditionMatrix(b11, b22, &MatrixTemp2);
      StrassenAlgorithm_omp(MatrixTemp1, MatrixTemp2, &p1);
    }
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      AdditionMatrix(a21, a22, &MatrixTemp1);
      StrassenAlgorithm_omp(MatrixTemp1, b11, &p2);
    }
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      SubtractionMatrix(b12, b22, &MatrixTemp1);
      StrassenAlgorithm_omp(a11, MatrixTemp1, &p3);
    }
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      SubtractionMatrix(b21, b11, &MatrixTemp1);
      StrassenAlgorithm_omp(a22, MatrixTemp1, &p4);
    }
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      AdditionMatrix(a11, a12, &MatrixTemp1);
      StrassenAlgorithm_omp(MatrixTemp1, b22, &p5);
    }
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      std::vector<double> MatrixTemp2(MatrixResultSize);
      SubtractionMatrix(a21, a11, &MatrixTemp1);
      AdditionMatrix(b11, b12, &MatrixTemp2);
      StrassenAlgorithm_omp(MatrixTemp1, MatrixTemp2, &p6);
    }
#pragma omp section
    {
      std::vector<double> MatrixTemp1(MatrixResultSize);
      std::vector<double> MatrixTemp2(MatrixResultSize);
      SubtractionMatrix(a12, a22, &MatrixTemp1);
      AdditionMatrix(b21, b22, &MatrixTemp2);
      StrassenAlgorithm_omp(MatrixTemp1, MatrixTemp2, &p7);
    }
  }

  std::vector<double> c11(MatrixResultSize);
  std::vector<double> c12(MatrixResultSize);
  std::vector<double> c21(MatrixResultSize);
  std::vector<double> c22(MatrixResultSize);

  AdditionMatrix(p3, p5, &c12);
  AdditionMatrix(p2, p4, &c21);
  std::vector<double> MatrixTemp1(MatrixResultSize);
  std::vector<double> MatrixTemp2(MatrixResultSize);
  AdditionMatrix(p1, p4, &MatrixTemp1);
  AdditionMatrix(MatrixTemp1, p7, &MatrixTemp2);
  SubtractionMatrix(MatrixTemp2, p5, &c11);

  std::vector<double> MatrixTemp3(MatrixResultSize);
  std::vector<double> MatrixTemp4(MatrixResultSize);
  AdditionMatrix(p1, p3, &MatrixTemp3);
  AdditionMatrix(MatrixTemp3, p6, &MatrixTemp4);
  SubtractionMatrix(MatrixTemp4, p2, &c22);

  RebuildMatrixOmp(MatrixResult, c11, c22, c12, c21, MatrixSize);
}
