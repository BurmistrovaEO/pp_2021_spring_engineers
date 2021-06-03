// Copyright 2021 Zlobin George
#ifndef MODULES_TASK_4_ZLOBIN_G_SPARSE_MATRIX_CCS_MULT_COMPLEX_STD_MATRIX_CCS_MULT_H_
#define MODULES_TASK_4_ZLOBIN_G_SPARSE_MATRIX_CCS_MULT_COMPLEX_STD_MATRIX_CCS_MULT_H_

#include <vector>
#include <complex>

class MatrixCCS {
 public:
    MatrixCCS(int nCollumns, int nRows, int nNotZero);
    MatrixCCS(int nCollumns, int nRows, std::vector<std::complex<int>> matrix);
    MatrixCCS(int nCollumns, int nRows,
              std::vector<int> collumnsIndexes,
              std::vector<int> rows,
              std::vector<std::complex<int>> values);

    void FillRandom(unsigned seed = 132, int min = -1000, int max = 1000);
    MatrixCCS Transpose() const;

    int GetNumCollumns() { return _nCollumns; }
    int GetNumRows() { return _nRows; }
    int GetNumNotZero() { return _nNotZero; }

    //  void Print();
 private:
    int _nCollumns;
    int _nRows;
    int _nNotZero;

    std::vector<int> _collumnsIndexes;
    std::vector<int> _rows;
    std::vector<std::complex<int>> _values;

    friend MatrixCCS SeqMultiplication(const MatrixCCS& A, const MatrixCCS& B);
    friend MatrixCCS STDMultiplication(const MatrixCCS& A, const MatrixCCS& B, int nThreads);
    friend bool operator==(const MatrixCCS& A, const MatrixCCS& B);
};

MatrixCCS SeqMultiplication(const MatrixCCS& A, const MatrixCCS& B);
MatrixCCS STDMultiplication(const MatrixCCS& A, const MatrixCCS& B, int nThreads = -1);
bool operator==(const MatrixCCS& A, const MatrixCCS& B);

#endif  // MODULES_TASK_4_ZLOBIN_G_SPARSE_MATRIX_CCS_MULT_COMPLEX_STD_MATRIX_CCS_MULT_H_
