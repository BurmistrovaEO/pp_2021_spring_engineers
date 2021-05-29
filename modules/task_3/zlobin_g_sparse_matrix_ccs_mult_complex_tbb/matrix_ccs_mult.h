// Copyright 2021 Zlobin George
#ifndef MODULES_TASK_3_ZLOBIN_G_SPARSE_MATRIX_CCS_MULT_COMPLEX_TBB_MATRIX_CCS_MULT_H_
#define MODULES_TASK_3_ZLOBIN_G_SPARSE_MATRIX_CCS_MULT_COMPLEX_TBB_MATRIX_CCS_MULT_H_

#include <tbb/tbb.h>
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
    friend MatrixCCS TBBMultiplication(const MatrixCCS& A, const MatrixCCS& B);
    friend bool operator==(const MatrixCCS& A, const MatrixCCS& B);

    friend class TBBMatrixMultiplicator;
};

class TBBMatrixMultiplicator {
 public:
    const MatrixCCS& AT;
    const MatrixCCS& B;

    std::vector<int>* const rows;
    std::vector<std::complex<int>>* const values;
    std::vector<int>* const collumnsIndexes;

    void operator()(const tbb::blocked_range<int> &range) const;

    TBBMatrixMultiplicator(const MatrixCCS& _AT,
                       const MatrixCCS& _B,
                       std::vector<int>* _rows,
                       std::vector<std::complex<int>>* _values,
                       std::vector<int>* _collumnsIndexes) :
    AT(_AT),
    B(_B),
    rows(_rows),
    values(_values),
    collumnsIndexes(_collumnsIndexes) {}
};

MatrixCCS SeqMultiplication(const MatrixCCS& A, const MatrixCCS& B);
MatrixCCS TBBMultiplication(const MatrixCCS& A, const MatrixCCS& B);
bool operator==(const MatrixCCS& A, const MatrixCCS& B);

#endif  // MODULES_TASK_3_ZLOBIN_G_SPARSE_MATRIX_CCS_MULT_COMPLEX_TBB_MATRIX_CCS_MULT_H_
