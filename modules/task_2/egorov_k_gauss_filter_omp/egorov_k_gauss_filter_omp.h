// Copyright 2020 Egorov Kirill
#ifndef MODULES_TASK_2_EGOROV_K_GAUSS_FILTER_OMP_EGOROV_K_GAUSS_FILTER_OMP_H_
#define MODULES_TASK_2_EGOROV_K_GAUSS_FILTER_OMP_EGOROV_K_GAUSS_FILTER_OMP_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int  sz);
int clamp(int value, int min, int max);
std::vector<float> calculateKernel3x3(float sigma);
int pColourCalculation(const std::vector<int> &img, int w, int h,
    int x, int y, std::vector<float> kernel);
std::vector<int> filterApply(const std::vector<int> &img,
    int w, int h, float sigma);
std::vector<int> filterApplyParallel(const std::vector<int> &img,
    int w, int h, float sigma);

#endif  // MODULES_TASK_2_EGOROV_K_GAUSS_FILTER_OMP_EGOROV_K_GAUSS_FILTER_OMP_H_
