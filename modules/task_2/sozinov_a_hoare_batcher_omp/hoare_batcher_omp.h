// Copyright 2021 Sozinov Alex

#ifndef MODULES_TASK_2_SOZINOV_A_HOARE_BATCHER_OMP_HOARE_BATCHER_OMP_H_
#define MODULES_TASK_2_SOZINOV_A_HOARE_BATCHER_OMP_HOARE_BATCHER_OMP_H_

#include <vector>

enum EvenOdd {
  Even,
  Odd
};

std::vector<double> GetRandomVector(int size);
void Sort(std::vector<double>* vector, int begin, int end);
void SeqSort(std::vector<double>* vector);
void EvenOddSplit(std::vector<double>* res, const std::vector<double>& left,
  const std::vector<double>& right, EvenOdd type);
void BatcherMerge(std::vector<double>* res, const std::vector<double>& left, const std::vector<double>& right);
void ParSort(std::vector<double>* vector, unsigned int numThreads = 1);
void FillOffset(std::vector<int>* offset, const int size, const int count);
unsigned int GetMaxNumThreads(const unsigned int numThread);

#endif  // MODULES_TASK_2_SOZINOV_A_HOARE_BATCHER_OMP_HOARE_BATCHER_OMP_H_
