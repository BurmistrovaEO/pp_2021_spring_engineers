// Copyright 2021 Ivanov Yaroslav
#ifndef MODULES_TASK_2_IVANOV_Y_DOUBLE_SORT_OPM_SORT_DOUBLE_H_
#define MODULES_TASK_2_IVANOV_Y_DOUBLE_SORT_OPM_SORT_DOUBLE_H_

#include <vector>

std::vector<double> BubbleForCompare(std::vector<double> inputNumbers, const int size);
std::vector<double> generateDouble(const int n);
std::vector<double> loc_sort(std::vector<double> inputNumers, const int size);
std::vector<double> transit(std::vector<double> inputNumers, std::vector<double> loc, const int size, int add);
std::vector<double> loc_sort_omp(std::vector<double> inputNumbers, const int size);
std::vector<double> transit_omp(std::vector<double> inputNumbers, std::vector<double> loc, const int size, int add);

#endif  // MODULES_TASK_2_IVANOV_Y_DOUBLE_SORT_OPM_SORT_DOUBLE_H_
