// Copyright 2021 Baulin Mikhail

#ifndef MODULES_TASK_2_BAULIN_M_RADIX_SORT_SIMPLE_MERGE_OMP_RADIX_SORT_H_
#define MODULES_TASK_2_BAULIN_M_RADIX_SORT_SIMPLE_MERGE_OMP_RADIX_SORT_H_
#include <vector>

std::vector<double> GetRandomVector(int size);
std::vector<double> RadixSort(const std::vector<double>& vec);
std::vector<double> Merge(const std::vector<double>& a,
                          const std::vector<double>& b);

#endif  // MODULES_TASK_2_BAULIN_M_RADIX_SORT_SIMPLE_MERGE_OMP_RADIX_SORT_H_
