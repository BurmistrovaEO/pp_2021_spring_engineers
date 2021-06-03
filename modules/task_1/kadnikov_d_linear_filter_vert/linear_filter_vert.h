// Copyright 2021 Kadnikov Denis
#ifndef MODULES_TASK_1_KADNIKOV_D_LINEAR_FILTER_VERT_LINEAR_FILTER_VERT_H_
#define MODULES_TASK_1_KADNIKOV_D_LINEAR_FILTER_VERT_LINEAR_FILTER_VERT_H_

#include <ctime>
#include <random>
#include <vector>

template <class T>
T Clamp(T value, T up, T low) {
  if (value < low) {
    return low;
  }
  if (value > up) {
    return up;
  }
  return value;
}

std::vector<double> gaussianFilterApply(const std::vector<double>& vec, int r,
                                        int c, int radius, int sigma);
std::vector<double> gaussianCalculate(int radius, int sigma);
std::vector<double> createRandomVector(int r, int c);
std::vector<double> transposition(const std::vector<double>& vec, int r, int c);

#endif  // MODULES_TASK_1_KADNIKOV_D_LINEAR_FILTER_VERT_LINEAR_FILTER_VERT_H_
