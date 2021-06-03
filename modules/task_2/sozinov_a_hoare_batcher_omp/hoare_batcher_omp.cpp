// Copyright 2021 Sozinov Alex

#include "../../../modules/task_2/sozinov_a_hoare_batcher_omp/hoare_batcher_omp.h"
#include <omp.h>
#include <algorithm>
#include <random>
#include <utility>
#include <iostream>
#include <vector>

std::vector<double> GetRandomVector(int size) {
  std::vector<double> genVec;
  std::random_device dev;
  std::mt19937 ger(dev());
  std::uniform_real_distribution<> realDist(-100, 100);

  for (int index = 0; index < size; ++index) {
    genVec.push_back(realDist(ger));
  }
  return genVec;
}

void Sort(std::vector<double>* vector, int begin, int end) {
  int first = begin, last = end;
  double mid = (*vector)[begin + (end - begin) / 2];

  do {
    while ((*vector)[first] < mid)
      ++first;
    while ((*vector)[last] > mid)
      --last;

    if (first <= last) {
      if (first < last)
        std::swap((*vector)[first], (*vector)[last]);
      ++first;
      --last;
    }
  } while (first <= last);

  if (first < end)
    Sort(vector, first, end);
  if (begin < last)
    Sort(vector, begin, last);
}

void SeqSort(std::vector<double>* vector) {
  std::vector<double> left(vector->begin(), vector->begin() + vector->size() / 2);
  std::vector<double> right(vector->begin() + vector->size() / 2, vector->end());
  Sort(&left, 0, left.size() - 1);
  Sort(&right, 0, right.size() - 1);

  vector->clear();
  BatcherMerge(vector, left, right);
}

void EvenOddSplit(std::vector<double>* res, const std::vector<double>& left,
  const std::vector<double>& right, EvenOdd type) {
  size_t leftIndex, rightIndex;
  if (type == EvenOdd::Even) {
    leftIndex = 0;
    rightIndex = 0;
  } else {
    leftIndex = 1;
    rightIndex = 1;
  }

  while ((rightIndex < right.size()) && (leftIndex < left.size())) {
    if (left[leftIndex] <= right[rightIndex]) {
      (*res).push_back(left[leftIndex]);
      leftIndex += 2;
    } else {
      (*res).push_back(right[rightIndex]);
      rightIndex += 2;
    }
  }

  if (leftIndex >= left.size()) {
    for (size_t i = rightIndex; i < right.size(); i += 2) {
      (*res).push_back(right[i]);
    }
  } else {
    for (size_t i = leftIndex; i < left.size(); i += 2) {
      (*res).push_back(left[i]);
    }
  }
}

void BatcherMerge(std::vector<double>* res, const std::vector<double>& left, const std::vector<double>& right) {
  std::vector<double> even, odd;
  size_t leftB = res->size();
  even.reserve((left.size() + 1) / 2 + (right.size() + 1) / 2);
  odd.reserve(left.size() / 2 + right.size() / 2);
  EvenOddSplit(&even, left, right, EvenOdd::Even);
  EvenOddSplit(&odd, left, right, EvenOdd::Odd);

  size_t indexL = 0, indexR = 0;
  while (indexL < even.size() && indexR < odd.size()) {
    if (even[indexL] < odd[indexR]) {
      (*res).push_back(even[indexL]);
      ++indexL;
    } else {
      (*res).push_back(odd[indexR]);
      ++indexR;
    }
  }

  if (indexR < odd.size())
    (*res).insert((*res).end(), odd.begin() + indexR, odd.end());
  else if (indexL < even.size())
    (*res).insert((*res).end(), even.begin() + indexL, even.end());

  for (size_t i = leftB; i < (*res).size() - 1; ++i) {
    if ((*res)[i] > (*res)[i + 1]) {
      std::swap((*res)[i], (*res)[i + 1]);
    }
  }
}

void ParSort(std::vector<double>* vector, unsigned int numThreads) {
  numThreads = GetMaxNumThreads(numThreads);
  if (vector->size() == 1)
    return;
  if (vector->size() < numThreads) {
    Sort(vector, 0, vector->size() - 1);
    return;
  }

  std::vector<int> offset(numThreads, 0);
  std::vector<double> resVec = std::vector<double>(vector->size(), 0), locVec;
  FillOffset(&offset, vector->size(), numThreads);

  omp_set_num_threads(numThreads);

#pragma omp parallel private(locVec) shared(vector, offset)
  {
    int numThread = omp_get_thread_num();
    locVec = std::vector<double>(offset[numThread + 1] - offset[numThread], 0);
    std::copy(vector->begin() + offset[numThread], vector->begin() + offset[numThread + 1], locVec.begin());
    Sort(&locVec, 0, static_cast<int>(locVec.size() - 1));

#pragma omp critical
    {
      std::copy(locVec.begin(), locVec.end(), resVec.begin() + offset[numThread]);
    }
  }

  unsigned int mergeNumThreads = (numThreads + numThreads % 2) / 2;
  if (mergeNumThreads < 1 || numThreads == 1) {
    *vector = std::vector<double>(resVec);
    return;
  }

  while (mergeNumThreads >= 1) {
    omp_set_num_threads(mergeNumThreads);
    std::vector<double> left, right;
    vector->clear();

#pragma omp parallel private(left, right) shared(resVec, vector)
    {
      int numThread = omp_get_thread_num();
      int leftB = offset[2 * numThread], rightB = offset[2 * numThread + 2];
      int mid = offset[2 * numThread + 1];

      std::copy(resVec.begin() + leftB, resVec.begin() + mid, std::back_inserter(left));
      std::copy(resVec.begin() + mid, resVec.begin() + rightB, std::back_inserter(right));
#pragma omp critical
      {
        BatcherMerge(vector, left, right);
      }
    }
    offset = std::vector<int>(mergeNumThreads);
    FillOffset(&offset, vector->size(), mergeNumThreads);
    resVec = std::vector<double>(*vector);
    mergeNumThreads /= 2;
  }
}

void FillOffset(std::vector<int>* offset, const int size, const int count) {
  int del = size % count;
  (*offset)[0] = 0;
  for (int i = 1; i < count; ++i) {
    int curSize = (i < del) ? size / count + 1 : size / count;
    (*offset)[i] = offset->at(i - 1) + curSize;
  }
  offset->push_back(size);
}

unsigned int GetMaxNumThreads(const unsigned int numThread) {
  unsigned int res = 1;
  while (res * 2 <= numThread)
    res *= 2;
  return res;
}
