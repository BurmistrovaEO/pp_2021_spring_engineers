// Copyright 2021 Pasukhin Dmitry

#include "../../../modules/task_2/pasukhin_d_crs/crs.h"

#include <complex>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

//=======================================================================
// function : TransformToNorm
// purpose  :
//=======================================================================

std::vector<std::complex<double>> TransformToNorm(const Matrix& A) {
  std::vector<std::complex<double>> Res(A.Lenght * A.Lenght, 0.);
  for (size_t i = 0; i < A.Lenght; ++i) {
    const size_t i1 = A.RowInd.at(i);
    const size_t i2 = A.RowInd.at(i + 1) - 1;
    for (size_t j = i1; j <= i2; j++) {
      Res.at(i * A.Lenght + A.Column.at(j)) = A.Values.at(j);
    }
  }
  return Res;
}

//=======================================================================
// function : GenerateCRS
// purpose  :
//=======================================================================

Matrix GenerateCRS(const size_t lenght) {
  Matrix Res;
  if (lenght < 1) {
    return Res;
  }

  std::mt19937 generator(static_cast<unsigned int>(time(0)));
  size_t ValueInRow = (lenght > 1) ? 1 + generator() % (lenght / 2) : 1;
  size_t VCount = ValueInRow * lenght;

  Res.Lenght = lenght;
  Res.VCount = VCount;
  Res.RowInd = std::vector<size_t>(lenght + 1);
  Res.Column = std::vector<size_t>(VCount);
  Res.Values = std::vector<std::complex<double>>(VCount);

  for (size_t i = 0; i < lenght; i++) {
    for (size_t j = 0; j < ValueInRow; j++) {
      bool isValue = false;
      do {
        isValue = false;
        Res.Column.at(i * ValueInRow + j) =
            (lenght > 1) ? 1 + generator() % (lenght - 1) : 0;
        for (size_t k = 0; k < j; k++) {
          if (Res.Column.at(i * ValueInRow + j) ==
              Res.Column.at(i * ValueInRow + k)) {
            isValue = true;
          }
        }
      } while (isValue);
    }
    for (size_t j = 0; j < ValueInRow - 1; j++) {
      for (size_t k = 0; k < ValueInRow - 1; k++) {
        if (Res.Column.at(i * ValueInRow + k) >
            Res.Column.at(i * ValueInRow + k + 1)) {
          std::swap(Res.Column.at(i * ValueInRow + k),
                    Res.Column.at(i * ValueInRow + k + 1));
        }
      }
    }
  }

  for (auto& iter : Res.Values) {
    iter = std::complex<double>(generator() % 110, generator() % 110);
  }

  size_t Count = 0;
  for (auto& iter : Res.RowInd) {
    iter = Count;
    Count += ValueInRow;
  }

  return Res;
}

//=======================================================================
// function : MultNorm
// purpose  :
//=======================================================================

std::vector<std::complex<double>> MultNorm(
    const std::vector<std::complex<double>>& A,
    const std::vector<std::complex<double>>& B, const size_t lenght) {
  std::vector<std::complex<double>> Res(A.size());
  for (size_t i = 0; i < lenght; i++) {
    for (size_t j = 0; j < lenght; j++) {
      std::complex<double> sum = 0;
      for (size_t k = 0; k < lenght; k++) {
        sum += A.at(i * lenght + k) * B.at(k * lenght + j);
      }
      Res.at(i * lenght + j) = sum;
    }
  }
  return Res;
}

//=======================================================================
// function : MultCRS
// purpose  :
//=======================================================================

Matrix MultCRS(const Matrix& A, const Matrix& b) {
  const size_t N = A.Lenght;
  const size_t N2 = b.Lenght;
  Matrix Res, B;

  B.Lenght = b.Lenght;
  B.VCount = b.VCount;
  B.Column = std::vector<size_t>(b.VCount);
  B.RowInd = std::vector<size_t>(b.Lenght + 1);
  B.Values = std::vector<std::complex<double>>(b.VCount);

  for (const auto& iter : b.Column) {
    B.RowInd.at(iter + 1)++;
  }

  size_t Index_tmp = 0;
  for (auto& iter : B.RowInd) {
    const size_t tmp = iter;
    iter = Index_tmp;
    Index_tmp = Index_tmp + tmp;
  }

  for (size_t i = 0; i < b.Lenght; i++) {
    size_t j1 = b.RowInd.at(i);
    size_t j2 = b.RowInd.at(i + 1);
    size_t Col = i;
    for (size_t j = j1; j < j2; j++) {
      std::complex<double> V = b.Values.at(j);
      size_t Row_index = b.Column.at(j);
      size_t IIndex = B.RowInd.at(Row_index + 1);
      B.Values.at(IIndex) = V;
      B.Column.at(IIndex) = Col;
      B.RowInd.at(Row_index + 1)++;
    }
  }

  Res.Lenght = N;
  Res.RowInd.push_back(0);
  for (size_t i = 0; i < N; i++) {
    size_t rowNZ = 0;
    for (size_t j = 0; j < N2; j++) {
      std::complex<double> sum = 0;
      for (size_t k = A.RowInd.at(i); k < A.RowInd.at(i + 1); k++) {
        for (size_t l = B.RowInd.at(j); l < B.RowInd.at(j + 1); l++) {
          if (A.Column.at(k) == B.Column.at(l)) {
            sum += A.Values.at(k) * B.Values.at(l);
            break;
          }
        }
      }

      if (abs(sum) > std::numeric_limits<double>::epsilon()) {
        Res.Column.push_back(j);
        Res.Values.push_back(sum);
        rowNZ++;
      }
    }
    Res.RowInd.push_back(rowNZ + Res.RowInd.at(i));
  }
  Res.VCount = Res.Column.size();
  return Res;
}

//=======================================================================
// function : CompareCRS
// purpose  :
//=======================================================================

bool Compare_Vector(const std::vector<size_t>& A,
                    const std::vector<size_t>& B) {
  for (size_t i = 0; i < A.size(); i++) {
    if (A.at(i) != B.at(i)) {
      return false;
    }
  }
  return true;
}

bool Compare_Vector(const std::vector<std::complex<double>>& A,
                    const std::vector<std::complex<double>>& B) {
  for (size_t i = 0; i < A.size(); i++) {
    if (abs(A.at(i) - B.at(i)) > 1e-8) {
      return false;
    }
  }
  return true;
}

bool CompareCRS(const Matrix& A, const Matrix& B) {
  if (A.Lenght != B.Lenght || A.VCount != B.VCount ||
      !Compare_Vector(A.Column, B.Column) ||
      !Compare_Vector(A.Values, B.Values) ||
      !Compare_Vector(A.RowInd, B.RowInd)) {
    return false;
  }
  return true;
}

//=======================================================================
// function : MultCRS_OMP
// purpose  :
//=======================================================================

Matrix MultCRS_OMP(const Matrix& A, const Matrix& b) {
  const size_t N = A.Lenght;
  const size_t N2 = b.Lenght;
  const int size = static_cast<int>(b.VCount);
  Matrix Res, B;

  B.Lenght = N2;
  B.VCount = N2;
  B.Column = std::vector<size_t>(b.VCount);
  B.RowInd = std::vector<size_t>(N2 + 2);
  B.Values = std::vector<std::complex<double>>(b.VCount);

  for (int i = 0; i < size; i++) {
    ++B.RowInd.at(b.Column.at(i) + 2);
  }

#pragma omp parallel for ordered
  for (int i = 2; i < static_cast<int>(N2) + 2; i++) {
#pragma omp ordered
    B.RowInd.at(i) += B.RowInd.at(i - 1);
  }

#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(N2); ++i) {
    const int end = static_cast<int>(b.RowInd.at(i + 1));
    for (int j = static_cast<int>(b.RowInd.at(i)); j < end; ++j) {
      {
        const size_t old_index = b.Column.at(j) + 1;
        size_t new_index = 0;
#pragma omp critical
        {
          new_index = B.RowInd.at(old_index);
          B.RowInd.at(old_index)++;
        }
        B.Values.at(new_index) = b.Values.at(j);
        B.Column.at(new_index) = i;
      }
    }
  }
  B.RowInd.pop_back();

  Res.Lenght = N;
  Res.RowInd = std::vector<size_t>(N2 + 1);
  std::vector<std::pair<size_t, std::complex<double>>> aColCalues(N * N);

#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(N); i++) {
    size_t rowNZ = 0;
    for (size_t j = 0; j < N2; j++) {
      std::complex<double> sum = 0;
      size_t end_k = A.RowInd.at(i + 1);
      size_t end_l = B.RowInd.at(j + 1);
      for (size_t k = A.RowInd.at(i); k < end_k; k++) {
        for (size_t l = B.RowInd.at(j); l < end_l; l++) {
          if (A.Column.at(k) == B.Column.at(l)) {
            sum += A.Values.at(k) * B.Values.at(l);
            break;
          }
        }
      }

      if (abs(sum) > 1e-8) {
        std::pair<size_t, std::complex<double>> aPair{j, sum};
        std::swap(aColCalues.at(i * N + rowNZ), aPair);
        rowNZ++;
      }
    }
    Res.RowInd.at(i + 1) = rowNZ;
  }
#pragma omp parallel for ordered
  for (int i = 1; i < static_cast<int>(N2) + 1; i++) {
#pragma omp ordered
    Res.RowInd.at(i) += Res.RowInd.at(i - 1);
  }
  const size_t res_size = Res.RowInd.back();
  Res.Column = std::vector<size_t>(res_size);
  Res.Values = std::vector<std::complex<double>>(res_size);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(N); i++) {
    const size_t start_ind = Res.RowInd.at(i);
    const size_t end_ind = Res.RowInd.at(i + 1);
    size_t own_ind = 0;
    for (size_t j = start_ind; j < end_ind; j++, own_ind++) {
      const std::pair<size_t, std::complex<double>> aPair(
          std::move(aColCalues.at(i * N + own_ind)));
      Res.Column.at(j) = aPair.first;
      Res.Values.at(j) = aPair.second;
    }
  }
  Res.VCount = Res.Column.size();
  return Res;
}
