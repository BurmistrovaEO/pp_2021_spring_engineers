// Copyright 2021 Igor Rukhovich
#ifndef MODULES_TASK_2_RUKHOVICH_I_CSR_MULT_DOUBLE_CSR_MULT_DOUBLE_H_
#define MODULES_TASK_2_RUKHOVICH_I_CSR_MULT_DOUBLE_CSR_MULT_DOUBLE_H_

#include <omp.h>

#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

class RandomDouble {
 public:
  static double Next() {
    static RandomDouble rand = RandomDouble();
    return rand.dist_(rand.gen_);
  }

 private:
  RandomDouble() : gen_(std::random_device()()), dist_(-1e3, 1e3) {
  }

  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> dist_;
};

template <class T>
bool Compare(const T &lhs, const T &rhs) {
  static float float_eps = 1e-3;
  static double double_eps = 1e-6;
  if (std::is_same<T, float>::value) {
    return (lhs + float_eps > rhs) && (lhs - float_eps < rhs);
  }
  if (std::is_same<T, double>::value) {
    return (lhs + double_eps > rhs) && (lhs - double_eps < rhs);
  }
  return lhs == rhs;
}

template <class ValueType, typename UIntType = uint16_t>
class CSRMatrixOMP {
 public:
  CSRMatrixOMP() = delete;

  CSRMatrixOMP(const CSRMatrixOMP &other) = default;
  CSRMatrixOMP(CSRMatrixOMP &&other) = default;

  explicit CSRMatrixOMP(UIntType height, UIntType width)
      : num_cols_(width), counts_(height + 1, 0) {
  }

  explicit CSRMatrixOMP(UIntType height, UIntType width, const std::vector<ValueType> &matrix)
      : num_cols_(width), counts_(height + 1) {
    if (static_cast<UIntType>(matrix.size()) != height * width) {
      throw std::runtime_error("Init matrix must consist of height*width elements");
    }
    counts_[0] = 0;
    for (UIntType row = 0; row < height; ++row) {
      counts_[row + 1] = counts_[row];
      for (UIntType col = 0; col < width; ++col) {
        auto cur_val = matrix[width * row + col];
        if (cur_val) {
          vals_.emplace_back(cur_val);
          cols_.emplace_back(col);
          ++counts_[row + 1];
        }
      }
    }
  }

  CSRMatrixOMP &operator=(CSRMatrixOMP other) {
    Swap(other);
    return *this;
  }

  void Swap(CSRMatrixOMP &other) {
    std::swap(vals_, other.vals_);
    std::swap(cols_, other.cols_);
    std::swap(counts_, other.counts_);
    std::swap(num_cols_, other.num_cols_);
  }

  bool operator==(const CSRMatrixOMP &rhs) const {
    const CSRMatrixOMP &lhs = *this;
    if (lhs.vals_.size() != rhs.vals_.size() || lhs.counts_.size() != rhs.counts_.size() ||
        lhs.num_cols_ != rhs.num_cols_) {
      return false;
    }
    for (UIntType i = 1; i < lhs.counts_.size(); ++i) {
      if (lhs.counts_[i] != rhs.counts_[i]) {
        return false;
      }
    }
    for (UIntType i = 0; i < lhs.vals_.size(); ++i) {
      if (lhs.cols_[i] != rhs.cols_[i] || !Compare(lhs.vals_[i], rhs.vals_[i])) {
        return false;
      }
    }
    return true;
  }

  CSRMatrixOMP &operator*=(const CSRMatrixOMP &other) {
    CSRMatrixOMP rhs = other.GetTransposed();
    UIntType res_width = other.num_cols_;
    UIntType res_height = counts_.size() - 1;
    std::vector<ValueType> res_mat(res_width * res_height);
#pragma omp parallel for schedule(dynamic)
    for (int32_t row = 0; row < static_cast<int32_t>(res_height); ++row) {
      for (UIntType col = 0; col < res_width; ++col) {
        UIntType l_row_cur = counts_[row], r_row_cur = rhs.counts_[col];
        ValueType cur_val = 0;
        while (l_row_cur < counts_[row + 1] && r_row_cur < rhs.counts_[col + 1]) {
          if (cols_[l_row_cur] < rhs.cols_[r_row_cur]) {
            ++l_row_cur;
          } else if (cols_[l_row_cur] > rhs.cols_[r_row_cur]) {
            ++r_row_cur;
          } else {
            cur_val += vals_[l_row_cur++] * rhs.vals_[r_row_cur++];
          }
        }
        res_mat[row * res_width + col] = cur_val;
      }
    }

    num_cols_ = res_width;
    vals_.resize(0);
    cols_.resize(0);
    counts_.resize(res_height + 1);
    counts_[0] = 0;
    for (UIntType row = 0; row < res_height; ++row) {
      counts_[row + 1] = counts_[row];
      for (UIntType col = 0; col < res_width; ++col) {
        auto cur_val = res_mat[res_width * row + col];
        if (cur_val) {
          vals_.emplace_back(cur_val);
          cols_.emplace_back(col);
          ++counts_[row + 1];
        }
      }
    }
    return *this;
  }

 protected:
  // We can think of it as converting CSR to CSC and back again
  const CSRMatrixOMP GetTransposed() const {
    CSRMatrixOMP other(num_cols_, counts_.size() - 1);
    for (UIntType i = 0; i < cols_.size(); ++i) {
      ++other.counts_[cols_[i]];
    }
    for (UIntType col = 0, sum = 0; col < other.counts_.size(); ++col) {
      UIntType tmp = other.counts_[col];
      other.counts_[col] = sum;
      sum += tmp;
    }
    other.cols_.resize(other.counts_.back());
    other.vals_.resize(other.counts_.back());

    for (UIntType row = 0; static_cast<UIntType>(row + 1) < counts_.size(); ++row) {
      for (UIntType cnt = counts_[row]; cnt < counts_[row + 1]; ++cnt) {
        UIntType col = cols_[cnt];
        UIntType dest_place = other.counts_[col];
        other.vals_[dest_place] = vals_[cnt];
        other.cols_[dest_place] = row;
        ++other.counts_[col];
      }
    }
    for (UIntType col = 0, last = 0; col <= num_cols_; ++col) {
      UIntType tmp = other.counts_[col];
      other.counts_[col] = last;
      last = tmp;
    }
    return other;
  }

  UIntType num_cols_;
  std::vector<ValueType> vals_;
  std::vector<UIntType> cols_;
  std::vector<UIntType> counts_;
};

#endif  // MODULES_TASK_2_RUKHOVICH_I_CSR_MULT_DOUBLE_CSR_MULT_DOUBLE_H_
