// Copyright 2021 Sozinov Alex

#include <algorithm>
#include <random>
#include <utility>
#include <iostream>
#include <vector>
#include "../../../3rdparty/unapproved/unapproved.h"
#include "../../../modules/task_4/sozinov_a_hoare_batcher_std/hoare_batcher_std.h"

std::vector<double> GetRandomVector(int size) {
  std::vector<double> genVec;
  std::random_device dev;
  std::mt19937 ger(dev());
  std::uniform_real_distribution<> realDist(-100, 100);

  genVec.reserve(size);
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

void ParSort(std::vector<double>* vector, unsigned int numThreads) {
  if (vector->size() == 1)
    return;
  if (vector->size() < numThreads) {
    Sort(vector, 0, vector->size() - 1);
    return;
  }

  std::vector<int> offset(numThreads, 0);
  std::vector<double> resVec = *vector, locVec;
  FillOffset(&offset, vector->size(), numThreads);

  std::vector<std::thread> threads;
  threads.reserve(numThreads);
  for (size_t index = 0; index < numThreads; ++index) {
      threads.emplace_back(std::thread([&offset, &resVec, index] {
          Sort(&resVec, offset[index], offset[index + 1] - 1);
      }));
  }
  for (size_t index = 0; index < threads.size(); ++index) {
      threads[index].join();
  }

  unsigned int mergeNumThreads = (numThreads + numThreads % 2) / 2;
  if (mergeNumThreads < 1 || numThreads == 1) {
    *vector = std::vector<double>(resVec);
    return;
  }

  vector->clear();
    for (unsigned int index = 0; index < numThreads; ++index) {
      std::vector<double> left, right;
      right = std::vector<double>(offset[index + 1] - offset[index], 0);
      std::copy(resVec.begin() + offset[index], resVec.begin() + offset[index + 1], right.begin());

      BatcherMergePar(vector, left, right);
    }
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
  even.reserve((left.size() + 1) / 2 + (right.size() + 1) / 2);
  odd.reserve(left.size() / 2 + right.size() / 2);
  EvenOddSplit(&even, left, right, EvenOdd::Even);
  EvenOddSplit(&odd, left, right, EvenOdd::Odd);

  size_t index = 0;
  for (; index < even.size() && index < odd.size(); ++index) {
    (*res).push_back(even[index]);
    (*res).push_back(odd[index]);
  }

  if (index < odd.size())
    (*res).insert((*res).end(), odd.begin() + index, odd.end());
  else if (index < even.size())
    (*res).insert((*res).end(), even.begin() + index, even.end());

  for (size_t i = 0; i < (*res).size() - 1; ++i) {
    if ((*res)[i] > (*res)[i + 1]) {
      std::swap((*res)[i], (*res)[i + 1]);
    }
  }
}

void BatcherMergePar(std::vector<double>* vector, const std::vector<double>& even, const std::vector<double>& odd) {
  size_t indexV = 0;
  std::vector<double> resVec;
  size_t indexL = 0, indexR = 0;

  while (indexL < vector->size() && indexR < odd.size()) {
    if (vector->at(indexL) < odd[indexR])
      resVec.push_back(vector->at(indexL++));
    else
      resVec.push_back(odd[indexR++]);
  }

  if (indexR < odd.size())
    resVec.insert(resVec.end(), odd.begin() + indexR, odd.end());
  if (indexV < vector->size())
    resVec.insert(resVec.end(), vector->begin() + indexL, vector->end());

  *vector = resVec;
  for (size_t i = 0; i < (*vector).size() - 1; ++i) {
    if ((*vector)[i] > (*vector)[i + 1]) {
      std::swap((*vector)[i], (*vector)[i + 1]);
    }
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
