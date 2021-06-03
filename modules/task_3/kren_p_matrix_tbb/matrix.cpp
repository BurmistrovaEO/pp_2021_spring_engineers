// Copyright 2021 Polina Kren

#include "../../modules/task_3/kren_p_matrix_tbb/matrix.h"

#include <math.h>
#include <tbb/task_group.h>
#include <tbb/tbb.h>

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

void RebuildMatrix(std::vector<double>* MatrixResult,
                   const std::vector<double>& A11,
                   const std::vector<double>& A22,
                   const std::vector<double>& A12,
                   const std::vector<double>& A21, const int MatrixSize) {
  const int MatrixResultSize = MatrixSize * 2;
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

void StrassenAlgorithm(const std::vector<double>& matrix1,
                       const std::vector<double>& matrix2,
                       std::vector<double>* result) {
  int MatrixSize = sqrt(result->size());
  if (MatrixSize <= 2) {
    MultiplicationMatrix(matrix1, matrix2, result, MatrixSize);
    return;
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

  std::vector<double> p1(MatrixResultSize);
  StrassenAlgorithm(AdditionMatrix(A11, A22, MatrixSize),
                    AdditionMatrix(B11, B22, MatrixSize), &p1);
  std::vector<double> p2(MatrixResultSize);
  StrassenAlgorithm(AdditionMatrix(A21, A22, MatrixSize), B11, &p2);
  std::vector<double> p3(MatrixResultSize);
  StrassenAlgorithm(A11, SubtractionMatrix(B12, B22, MatrixSize), &p3);
  std::vector<double> p4(MatrixResultSize);
  StrassenAlgorithm(A22, SubtractionMatrix(B21, B11, MatrixSize), &p4);
  std::vector<double> p5(MatrixResultSize);
  StrassenAlgorithm(AdditionMatrix(A11, A12, MatrixSize), B22, &p5);
  std::vector<double> p6(MatrixResultSize);
  StrassenAlgorithm(SubtractionMatrix(A21, A11, MatrixSize),
                    AdditionMatrix(B11, B12, MatrixSize), &p6);
  std::vector<double> p7(MatrixResultSize);
  StrassenAlgorithm(SubtractionMatrix(A12, A22, MatrixSize),
                    AdditionMatrix(B21, B22, MatrixSize), &p7);

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

  RebuildMatrix(result, c11, c22, c12, c21, MatrixSize);
}

void StrassenAlgorithm_tbb(const std::vector<double>& matrix1,
                           const std::vector<double>& matrix2,
                           std::vector<double>* matrix3) {
  int MatrixSize = sqrt(matrix1.size());
  if (MatrixSize <= 2) {
    MultiplicationMatrix(matrix1, matrix2, matrix3, MatrixSize);
    return;
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

  std::vector<double> MainP1(MatrixResultSize);
  std::vector<double>* p1 = &MainP1;
  std::vector<double> MainP2(MatrixResultSize);
  std::vector<double>* p2 = &MainP2;
  std::vector<double> MainP3(MatrixResultSize);
  std::vector<double>* p3 = &MainP3;
  std::vector<double> MainP4(MatrixResultSize);
  std::vector<double>* p4 = &MainP4;
  std::vector<double> MainP5(MatrixResultSize);
  std::vector<double>* p5 = &MainP5;
  std::vector<double> MainP6(MatrixResultSize);
  std::vector<double>* p6 = &MainP6;
  std::vector<double> MainP7(MatrixResultSize);
  std::vector<double>* p7 = &MainP7;

  tbb::task_scheduler_init init(8);
  tbb::task_group tasks;
  tasks.run([=] {
    StrassenAlgorithm(AdditionMatrix(A11, A22, MatrixSize),
                      AdditionMatrix(B11, B22, MatrixSize), p1);
  });

  tasks.run([=] {
    StrassenAlgorithm(AdditionMatrix(A21, A22, MatrixSize), B11, p2);
  });

  tasks.run([=] {
    StrassenAlgorithm(A11, SubtractionMatrix(B12, B22, MatrixSize), p3);
  });

  tasks.run([=] {
    StrassenAlgorithm(A22, SubtractionMatrix(B21, B11, MatrixSize), p4);
  });

  tasks.run([=] {
    StrassenAlgorithm(AdditionMatrix(A11, A12, MatrixSize), B22, p5);
  });

  tasks.run([=] {
    StrassenAlgorithm(SubtractionMatrix(A21, A11, MatrixSize),
                      AdditionMatrix(B11, B12, MatrixSize), p6);
  });

  tasks.run([=] {
    StrassenAlgorithm(SubtractionMatrix(A12, A22, MatrixSize),
                      AdditionMatrix(B21, B22, MatrixSize), p7);
  });
  tasks.wait();

  std::vector<double> c11(MatrixResultSize);
  std::vector<double> c12(MatrixResultSize);
  std::vector<double> c21(MatrixResultSize);
  std::vector<double> c22(MatrixResultSize);

  AdditionMatrix(MainP3, MainP5, &c12);
  AdditionMatrix(MainP2, MainP4, &c21);
  AdditionMatrix(AdditionMatrix(MainP1, MainP4, MatrixSize),
                 SubtractionMatrix(MainP7, MainP5, MatrixSize), &c11);
  AdditionMatrix(SubtractionMatrix(MainP1, MainP2, MatrixSize),
                 AdditionMatrix(MainP3, MainP6, MatrixSize), &c22);

  RebuildMatrix(matrix3, c11, c22, c12, c21, MatrixSize);
}
