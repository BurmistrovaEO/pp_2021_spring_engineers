// Copyright 2021 Bulychev Andrey
#ifndef MODULES_TASK_3_BULYCHEV_A_SHELL_SORT_SIMPLE_TBB_SHELL_SORT_SIMPLE_TBB_H_
#define MODULES_TASK_3_BULYCHEV_A_SHELL_SORT_SIMPLE_TBB_SHELL_SORT_SIMPLE_TBB_H_

#include <vector>


std::vector<double> getRandomVector(int size);
std::vector<double> Shell_sort(const std::vector<double>& vec);
std::vector<double> Merge(const std::vector<double>& vec1, const std::vector<double>& vec2);
std::vector<double> Shell_Sort_tbb(const std::vector<double>& vec, int num_threads = 4);


#endif  // MODULES_TASK_3_BULYCHEV_A_SHELL_SORT_SIMPLE_TBB_SHELL_SORT_SIMPLE_TBB_H_
