// Copyright 2021 Zlobin George
#include "../../../modules/task_3/zlobin_g_sparse_matrix_ccs_mult_complex_tbb/matrix_ccs_mult.h"
#include <random>
#include <algorithm>

#include <iostream>

MatrixCCS::MatrixCCS(int nCollumns, int nRows, int nNotZero) {
    if (nCollumns <= 0 || nRows <= 0) {
        throw "Incorrect matrix size";
    } else if (nNotZero < 0 || nNotZero > nCollumns * nRows) {
        throw "Incorrect number of not zero elements";
    }

    _nCollumns = nCollumns;
    _nRows = nRows;
    _nNotZero = nNotZero;

    _collumnsIndexes.reserve(_nCollumns + 1);
    _rows.reserve(_nRows);
    _values.reserve(_nNotZero);
}

MatrixCCS::MatrixCCS(int nCollumns, int nRows,
                     std::vector<std::complex<int>> matrix) {
    if (nCollumns <= 0 || nRows <= 0) {
        throw "Incorrect matrix size";
    } else if (static_cast<int>(matrix.size()) != nCollumns * nRows) {
        throw "Size of matrix doesn't match given parameters";
    }

    _nCollumns = nCollumns;
    _nRows = nRows;

    _collumnsIndexes.reserve(_nCollumns + 1);
    _rows.reserve(_nRows);

    _collumnsIndexes.push_back(0);
    for (int j = 0; j < _nCollumns; j++) {
        for (int i = 0; i < _nRows; i++) {
            if (matrix[i * _nCollumns + j] != std::complex<int>(0, 0)) {
                _values.push_back(matrix[i * _nCollumns + j]);
                _rows.push_back(i);
            }
        }
        _collumnsIndexes.push_back(static_cast<int>(_values.size()));
    }

    _nNotZero = _collumnsIndexes[_nCollumns];
}

MatrixCCS::MatrixCCS(int nCollumns, int nRows,
                     std::vector<int> collumnsIndexes,
                     std::vector<int> rows,
                     std::vector<std::complex<int>> values) {
    if (nCollumns <= 0 || nRows <= 0) {
        throw "Incorrect matrix size";
    } else if (static_cast<int>(values.size()) > nCollumns * nRows) {
        throw "More elements when possible";
    } else if (static_cast<int>(collumnsIndexes.size()) != nCollumns + 1) {
        throw "Incorrect collumns indexes";
    } else if (std::find_if(rows.begin(), rows.end(),
               [&nRows](int i) { return i >= nRows; }) != rows.end()) {
        throw "Incorrect rows";
    }

    _nCollumns = nCollumns;
    _nRows = nRows;
    _nNotZero = static_cast<int>(values.size());
    _collumnsIndexes = collumnsIndexes;
    _rows = rows;
    _values = values;
}

MatrixCCS MatrixCCS::Transpose() const {
    MatrixCCS transposed{_nRows, _nCollumns, _nNotZero};

    transposed._collumnsIndexes.resize(_nRows + 1);
    for (int i = 0; i < _nNotZero; i++) {
        transposed._collumnsIndexes[_rows[i] + 1]++;
    }

    int tmp, sum = 0;
    for (int i = 1; i <= _nRows; i++) {
        tmp = transposed._collumnsIndexes[i];
        transposed._collumnsIndexes[i] = sum;
        sum += tmp;
    }

    transposed._rows.resize(_nNotZero);
    transposed._values.resize(_nNotZero);

    int j1, j2;
    int row, column, place;
    std::complex<int> value;
    for (int i = 0; i < _nCollumns; i++) {
        j1 = _collumnsIndexes[i];
        j2 = _collumnsIndexes[i + 1];
        row = i;
        for (int j = j1; j < j2; j++) {
            value = _values[j];
            column = _rows[j];
            place = transposed._collumnsIndexes[column + 1];
            transposed._rows[place] = row;
            transposed._values[place] = value;
            transposed._collumnsIndexes[column + 1]++;
        }
    }

    return transposed;
}

MatrixCCS SeqMultiplication(const MatrixCCS& A, const MatrixCCS& B) {
    if (A._nCollumns != B._nRows) {
        throw "Can't multiply matrix which not corresponding";
    }

    MatrixCCS AT = A.Transpose();

    std::vector<int>* rows = new std::vector<int>[B._nCollumns];
    std::vector<std::complex<int>>* values = new std::vector<std::complex<int>>[B._nCollumns];

    std::vector<int> collumnsIndexes(B._nCollumns + 1, 0);

    std::vector<int> tempRows(B._nRows);
    for (int i = 0; i < B._nCollumns; i++) {
        // Clear tempRows
        tempRows = std::vector<int>(B._nRows, -1);

        // Fill tempRows with rows indexes
        int ks = B._collumnsIndexes[i];
        int kf = B._collumnsIndexes[i + 1];
        for (int k = ks; k < kf; k++) {
            int row = B._rows[k];
            tempRows[row] = k;
        }

        // Multiply on every collumn of AT
        for (int j = 0; j < AT._nCollumns; j++) {
            std::complex<int> value{0, 0};

            int ls = AT._collumnsIndexes[j];
            int lf = AT._collumnsIndexes[j + 1];
            for (int l = ls; l < lf; l++) {
                int rowAT = AT._rows[l];
                int rowB = tempRows[rowAT];
                if (rowB != -1) {
                    value += B._values[rowB] * AT._values[l];
                }
            }

            if (value != std::complex<int>(0, 0)) {
                values[i].push_back(value);
                rows[i].push_back(j);
                collumnsIndexes[i]++;
            }
        }
    }

    int nNotZero = 0;
    for (int i = 0; i < B._nCollumns; i++) {
        int tmp = collumnsIndexes[i];
        collumnsIndexes[i] = nNotZero;
        nNotZero += tmp;
    }
    collumnsIndexes[B._nCollumns] = nNotZero;

    MatrixCCS product{B._nCollumns, A._nRows, nNotZero};
    product._collumnsIndexes = collumnsIndexes;
    for (int i = 0; i < B._nCollumns; i++) {
        product._rows.insert(product._rows.end(), rows[i].begin(), rows[i].end());
        product._values.insert(product._values.end(), values[i].begin(), values[i].end());
    }

    delete[] rows;
    delete[] values;

    return product;
}

MatrixCCS TBBMultiplication(const MatrixCCS& A, const MatrixCCS& B) {
    if (A._nCollumns != B._nRows) {
        throw "Can't multiply matrix which not corresponding";
    }

    MatrixCCS AT = A.Transpose();

    std::vector<int>* rows = new std::vector<int>[B._nCollumns];
    std::vector<std::complex<int>>* values = new std::vector<std::complex<int>>[B._nCollumns];

    std::vector<int> collumnsIndexes(B._nCollumns + 1, 0);

    TBBMatrixMultiplicator functor(AT, B,
                          rows,
                          values,
                          &collumnsIndexes);
    tbb::task_scheduler_init init;
    tbb::parallel_for(tbb::blocked_range<int>(0, B._nCollumns, 100), functor);
    init.terminate();

    int nNotZero = 0;
    for (int i = 0; i < B._nCollumns; i++) {
        int tmp = collumnsIndexes[i];
        collumnsIndexes[i] = nNotZero;
        nNotZero += tmp;
    }
    collumnsIndexes[B._nCollumns] = nNotZero;

    MatrixCCS product{B._nCollumns, A._nRows, nNotZero};
    product._collumnsIndexes = collumnsIndexes;
    for (int i = 0; i < B._nCollumns; i++) {
        product._rows.insert(product._rows.end(), rows[i].begin(), rows[i].end());
        product._values.insert(product._values.end(), values[i].begin(), values[i].end());
    }

    delete[] rows;
    delete[] values;

    return product;
}

void TBBMatrixMultiplicator::operator()(const tbb::blocked_range<int> &range) const {
    for (int i = range.begin(); i < range.end(); i++) {
        // Clear tempRows
        std::vector<int> tempRows = std::vector<int>(B._nRows, -1);

        // Fill tempRows with rows indexes
        int ks = B._collumnsIndexes[i];
        int kf = B._collumnsIndexes[i + 1];
        for (int k = ks; k < kf; k++) {
            int row = B._rows[k];
            tempRows[row] = k;
        }

        // Multiply on every collumn of AT
        for (int j = 0; j < AT._nCollumns; j++) {
            std::complex<int> value{0, 0};

            int ls = AT._collumnsIndexes[j];
            int lf = AT._collumnsIndexes[j + 1];
            for (int l = ls; l < lf; l++) {
                int rowAT = AT._rows[l];
                int rowB = tempRows[rowAT];
                if (rowB != -1) {
                    value += B._values[rowB] * AT._values[l];
                }
            }

            if (value != std::complex<int>(0, 0)) {
                values[i].push_back(value);
                rows[i].push_back(j);
                (*collumnsIndexes)[i]++;
            }
        }
    }
}

void MatrixCCS::FillRandom(unsigned seed, int min, int max) {
    _collumnsIndexes.clear();
    _rows.clear();
    _values.clear();

    std::mt19937 random;
    random.seed(seed);

    std::complex<int> value;
    const double probBarrier = 1.0 * _nNotZero / (_nCollumns * _nRows);

    std::uniform_real_distribution<> probDistrib(0.0, 1.0);
    std::uniform_int_distribution<int> valueDistrib(min, max);

    _collumnsIndexes.push_back(0);
    for (int collumn = 0; collumn < _nCollumns; collumn++) {
        for (int row = 0; row < _nRows; row++) {
            if (static_cast<int>(_values.size()) == _nNotZero) {
                break;
            } else if (_nNotZero - static_cast<int>(_values.size()) + collumn * _nRows + row ==
                       _nCollumns * _nRows ||
                       (static_cast<int>(_values.size()) < _nNotZero &&
                        probDistrib(random) < probBarrier)) {
                do {
                    value = std::complex<int>(valueDistrib(random), valueDistrib(random));
                } while (value == std::complex<int>(0, 0));
                _values.push_back(value);
                _rows.push_back(row);
            }
        }
        _collumnsIndexes.push_back(static_cast<int>(_values.size()));
        if (static_cast<int>(_values.size()) == _nNotZero) {
            while (static_cast<int>(_collumnsIndexes.size()) < _nCollumns + 1) {
                _collumnsIndexes.push_back(0);
            }
            break;
        }
    }
}

// void MatrixCCS::Print() {
//     for (int i = 0; i < _nCollumns; i++) {
//         for (int k = _collumnsIndexes[i];
//                 k < _collumnsIndexes[i + 1]; k++) {
//             std::cout << "Column " << i << "\tRow " << _rows[k] <<
//                 "\tValueR" << _values[k].real() << "\tValueI" << _values[k].imag() << "\n";
//         }
//     }
// }

bool operator==(const MatrixCCS& A, const MatrixCCS& B) {
    return A._nCollumns == B._nCollumns &&
           A._nRows == B._nRows &&
           A._nNotZero == B._nNotZero &&
           A._collumnsIndexes == B._collumnsIndexes &&
           A._rows == B._rows &&
           A._values == B._values;
}
