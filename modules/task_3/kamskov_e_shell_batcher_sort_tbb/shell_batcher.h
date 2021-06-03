// Copyright 2021 Kamskov Eugene
#ifndef MODULES_TASK_3_KAMSKOV_E_SHELL_BATCHER_SORT_TBB_SHELL_BATCHER_H_
#define MODULES_TASK_3_KAMSKOV_E_SHELL_BATCHER_SORT_TBB_SHELL_BATCHER_H_

#include <tbb/tbb.h>
#include <vector>
#include <ctime>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <random>

std::vector<int> getRandomVector(int);
std::vector<std::vector<int>> splitV(const std::vector<int>&, size_t);
std::vector<int> shellSort(const std::vector<int>&, int);
std::vector<int> evenBatch(const std::vector<int>&, const std::vector<int>&);
std::vector<int> oddBatch(const std::vector<int>&, const std::vector<int>&);
std::vector<int> batchShellSort(const std::vector<int>&, int);
std::vector<int> batchShellSortTbb(const std::vector<int>&, const int, int s);
std::vector<int> mergeBatch(const std::vector<int>&, const std::vector<int>&);
void print(std::vector<int> const&);
#endif  // MODULES_TASK_3_KAMSKOV_E_SHELL_BATCHER_SORT_TBB_SHELL_BATCHER_H_
