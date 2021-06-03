// Copyright 2021 Vizgalov Anton

#ifndef MODULES_TASK_3_VIZGALOV_A_HOARE_SORT_SIMPLE_MERGE_HOARE_SORT_SIMPLE_MERGE_H_
#define MODULES_TASK_3_VIZGALOV_A_HOARE_SORT_SIMPLE_MERGE_HOARE_SORT_SIMPLE_MERGE_H_

#include <vector>

int partition(std::vector<int>* vec, int low, int high);
void hoareSort(std::vector<int>* vec, int low, int high);
std::vector<int> merge(const std::vector<int>& vec1, const std::vector<int>& vec2);

void runHoareSort(std::vector<int>* vec, int numThreads = 4);

std::vector<int> getRandomVector(int size, int type = 0);

#endif  // MODULES_TASK_3_VIZGALOV_A_HOARE_SORT_SIMPLE_MERGE_HOARE_SORT_SIMPLE_MERGE_H_
