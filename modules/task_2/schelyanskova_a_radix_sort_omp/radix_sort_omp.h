// Copyright 2021 Schelyanskova Anastasiia
#ifndef MODULES_TASK_2_SCHELYANSKOVA_A_RADIX_SORT_OMP_RADIX_SORT_OMP_H_
#define MODULES_TASK_2_SCHELYANSKOVA_A_RADIX_SORT_OMP_RADIX_SORT_OMP_H_

#include <omp.h>

#include <string>
#include <vector>

std::vector<int> getRandomVector(int size);
int get_max_power(std::vector<int> vect);
bool vector_sort(std::vector<int> vect);
std::vector<int> RadixSort(std::vector<int> vect, int size);
std::vector<int> RadixSortParallel(std::vector<int> vect, int size,
                                   int threads);

// std::vector<int> Merge(const std::vector<int> &vect_left,
//                      const std::vector<int> &vect_right);


#endif  // MODULES_TASK_2_SCHELYANSKOVA_A_RADIX_SORT_OMP_RADIX_SORT_OMP_H_
