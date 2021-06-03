// Copyright 2021 Kadnikov Denis
#include "../../../modules/task_1/kadnikov_d_linear_filter_vert/linear_filter_vert.h"

#include <algorithm>
#include <cmath>

std::vector<double> transposition(const std::vector<double>& vec, int r,
                                  int c) {
  std::vector<double> res(r * c);

  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      res[i + j * r] = vec[i * c + j];
    }
  }

  return res;
}

std::vector<double> createRandomVector(int r, int c) {
  if ((r <= 0) || (c <= 0)) throw 1;

  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));

  std::vector<double> res(r * c);

  for (int i = 0; i < r * c; i++) res[i] = gen() % 256;

  return res;
}

std::vector<double> gaussianCalculate(int rad, int sigma) {
  if ((rad <= 0) || (sigma <= 0)) throw 1;

  const int size = 2 * rad + 1;
  double normal = 0;

  std::vector<double> filter(size * size);

  for (int i = -rad; i <= rad; i++) {
    for (int j = -rad; j <= rad; j++) {
      int index = (i + rad) * size + j + rad;
      filter[index] = exp(-(i * i + j * j) / (sigma * sigma));
      normal += filter[index];
    }
  }

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) filter[i * size + j] /= normal;
  return filter;
}

std::vector<double> gaussianFilterApply(const std::vector<double>& vec, int r,
                                        int c, int rad, int sigma) {
  const unsigned int size = 2 * rad + 1;

  std::vector<double> resVec(r * c);
  std::vector<double> filter = gaussianCalculate(rad, sigma);

  for (int x = 0; x < r; x++) {
    for (int y = 0; y < c; y++) {
      double res = 0;

      for (int i = -rad; i <= rad; i++) {
        for (int j = -rad; j <= rad; j++) {
          int index = (i + rad) * size + j + rad;
          int x_ = Clamp(x + j, r - 1, 0);
          int y_ = Clamp(y + i, c - 1, 0);
          double value = vec[x_ * c + y_];
          res += value * filter[index];
        }
      }
      res = Clamp(static_cast<int>(res), 255, 0);
      resVec[x * c + y] = res;
    }
  }
  return resVec;
}
