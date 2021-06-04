// Copyright 2021 Solomakhin Sergey
#ifndef MODULES_TASK_3_SOLOMAKHIN_S_HOAR_SORT_SIMPLE_MERGE_HOAR_SORT_SIMPLE_MERGE_H_
#define MODULES_TASK_3_SOLOMAKHIN_S_HOAR_SORT_SIMPLE_MERGE_HOAR_SORT_SIMPLE_MERGE_H_

#include <vector>

std::vector<int> random_gen(int size, int type = 0);
void hoar_sort(std::vector<int>* arr, int left, int right);
int partition(std::vector<int>* arr, int left, int right);
bool checker(std::vector <int>* vec);
std::vector<int> merge(const std::vector<int>& a, const std::vector<int>& b);
void hoar_sort_tbb(std::vector<int>* vec, int numThreads = 4);

#endif  // MODULES_TASK_3_SOLOMAKHIN_S_HOAR_SORT_SIMPLE_MERGE_HOAR_SORT_SIMPLE_MERGE_H_
