// Copyright 2021 Schelyanskova Anastasiia
#ifndef MODULES_TASK_3_SCHELYANSKOVA_A_RADIX_SORT_TBB_BITWISE_SORT_SIMPLE_TBB_H_
#define MODULES_TASK_3_SCHELYANSKOVA_A_RADIX_SORT_TBB_BITWISE_SORT_SIMPLE_TBB_H_

#include <vector>


std::vector<int> getRandomVector(int size);
std::vector<int> RadixSort(std::vector<int> vect, int size);
std::vector<int> RadixSort_tbb(const std::vector<int>& v);
std::vector<int> Merge(const std::vector<int>& vect_left,
                          const std::vector<int>& vect_right);


#endif  // MODULES_TASK_3_SCHELYANSKOVA_A_RADIX_SORT_TBB_BITWISE_SORT_SIMPLE_TBB_H_
