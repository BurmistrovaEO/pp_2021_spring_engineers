// Copyright 2021 Kamskov Eugene
#ifndef MODULES_TASK_2_KAMSKOV_E_SHELL_BATCHER_SORT_OMP_SHELL_BATCHER_H_
#define MODULES_TASK_2_KAMSKOV_E_SHELL_BATCHER_SORT_OMP_SHELL_BATCHER_H_

#include <omp.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <utility>
#include <random>
#include <ctime>

std::vector<int> getRandomVector(int);
std::vector<int> shellSort(const std::vector<int>&, int size);
std::vector<int> evenBatch(const std::vector<int>&, const std::vector<int>&);
std::vector<int> oddBatch(const std::vector<int>&, const std::vector<int>&);
std::vector<int> mergeBatch(const std::vector<int>&, const std::vector<int>&);
std::vector<int> batchShellSortOmp(const std::vector<int>&, const int, int);
std::vector<int> mergeOmp(const std::vector<std::vector<int>>&, const int, int);
std::vector<int> batchShellSort(const std::vector<int>&, int);
void print(std::vector<int> const&);
#endif  // MODULES_TASK_2_KAMSKOV_E_SHELL_BATCHER_SORT_OMP_SHELL_BATCHER_H_
