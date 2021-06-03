// Copyright 2021 Evseev Alexander
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <stdexcept>
#include <vector>

#include "../../../modules/task_2/evseev_a_gauss_omp/evseev_a_gauss_omp.h"

int clamp(int value, int max, int min) {
  if (value < min) {
    return min;
  }
  if (value > max) {
    return max;
  }
  return value;
}
Matrix RandMatrix(int height, int width) {
  if ((height <= 0) || (width <= 0)) throw std::invalid_argument("Invalid size");
  std::mt19937 gen;
  gen.seed(static_cast<int>(time(0)));
  Matrix matrix(height, std::vector<double>(width));
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      matrix[i][j] = gen() % 100;
    }
  }
  return matrix;
}

Matrix GaussKernel(int R, double sigma) {
  const int size = 2 * R + 1;
  double norm = 0;
  Matrix kernel(size, std::vector<double>(size));
  for (int64_t i = -R; i <= R; i++) {
    for (int64_t j = -R; j <= R; j++) {
      kernel[i + R][j + R] =
          static_cast<double>(exp(-(i * i + j * j) / (sigma * sigma)));
      norm = norm + kernel[i + R][j + R];
    }
  }
  for (int64_t i = 0; i < size; i++) {
    for (int64_t j = 0; j < size; j++) {
      kernel[i][j] = kernel[i][j] / norm;
    }
  }
  return kernel;
}

Matrix SeqGauss(const Matrix& mat, int height, int width, int R, double sigma) {
  Matrix finMat(height, std::vector<double>(width));
  Matrix kernel = GaussKernel(R, sigma);
  for (int64_t x = 0; x < height; x++) {
    for (int64_t y = 0; y < width; y++) {
      int finValue = 0;
      for (int64_t i = -R; i <= R; i++) {
        for (int64_t j = -R; j <= R; j++) {
          double value = mat[x][y];
          finValue += value * kernel[i + R][j + R];
        }
      }
      finValue = clamp(finValue, 255, 0);
      finMat[x][y] = finValue;
    }
  }
  return finMat;
}

Matrix OmpGauss(const Matrix& mat, int height, int width, int R, double sigma) {
  Matrix finMat(height, std::vector<double>(width));
  Matrix kernel = GaussKernel(R, sigma);
#pragma omp parallel
  {
#pragma omp for collapse(2) schedule(static)
    for (int64_t x = 0; x < height; x++) {
      for (int64_t y = 0; y < width; y++) {
        int finValue = 0;
        for (int64_t i = -R; i <= R; i++) {
          for (int64_t j = -R; j <= R; j++) {
            double value = mat[x][y];
            finValue += value * kernel[i + R][j + R];
          }
        }
        finValue = clamp(finValue, 255, 0);
        finMat[x][y] = finValue;
      }
    }
  }
  return finMat;
}
