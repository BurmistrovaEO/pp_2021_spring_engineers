// Copyright 2021 Frolova Olga
#ifndef MODULES_TASK_2_FROLOVA_O_RADIX_SORT_BATCHER_MERGE_OMP_UPDATE_RADIX_SORT_BATCHER_MERGE_UPDATE_H_
#define MODULES_TASK_2_FROLOVA_O_RADIX_SORT_BATCHER_MERGE_OMP_UPDATE_RADIX_SORT_BATCHER_MERGE_UPDATE_H_

#include <vector>

std::vector<double> getRandomVector(int size);
std::vector<double> radixSortPos(std::vector<double> vec, int radix);
std::vector<double> checkVector(std::vector<double> vec);
double maxVal(std::vector<double> vec);
int countRadix(double num);
int countNegRadix(double num);
int getRemainder(double num);
std::vector<double> radixSortNeg(std::vector<double> vec, int radix);
std::vector<double> radixSort(std::vector<double> vec);
std::vector<double> radix_sort_batcher_omp(std::vector<double> vec,
                                           int num_threads);
void makeNetwork(int size);
void net(std::vector<int> vec);
void oddEvenMerge(std::vector<int> left, std::vector<int> right);

#endif  // MODULES_TASK_2_FROLOVA_O_RADIX_SORT_BATCHER_MERGE_OMP_UPDATE_RADIX_SORT_BATCHER_MERGE_UPDATE_H_
