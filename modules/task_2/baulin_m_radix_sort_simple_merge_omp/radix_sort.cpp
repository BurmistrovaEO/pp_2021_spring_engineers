// Copyright 2021 Baulin Mikhail

#include "../../../modules/task_2/baulin_m_radix_sort_simple_merge_omp/radix_sort.h"

#include <omp.h>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>

std::vector<double> GetRandomVector(int sz) {
  std::random_device dev;
  std::mt19937 gen(dev());
  std::vector<double> vec(sz);
  for (int i = 0; i < sz; i++) {
    vec[i] = gen() % 100;
  }
  return vec;
}

std::vector<double> CountingSort(std::vector<double> source,
                                 std::vector<double> tmp, const int offset) {
  auto const size = static_cast<int>(source.size());
  auto byteArr = (unsigned char*)source.data();
  auto acc = 0;
  int counters[256] = {0};

  for (auto i = 0; i < size; i++) {
    counters[byteArr[offset + 8 * i]]++;
  }

  for (auto i = 0; i < 256; i++) {
    int tmp = counters[i];
    counters[i] = acc;
    acc += tmp;
  }

  for (auto i = 0; i < size; i++) {
    auto byteIndex = offset + 8 * i;
    auto byte = byteArr[byteIndex];
    auto counter = counters[byte];
    tmp[counter] = source[i];
    counters[byte]++;
  }

  return tmp;
}

std::vector<double> SignSort(std::vector<double> source) {
  auto const offset = 7;
  auto const size = static_cast<int>(source.size());
  auto byteArr = (unsigned char*)source.data();
  auto acc = 0;

  auto result = std::vector<double>(source);
  int counters[256] = {0};

  for (auto i = 0; i < size; i++) {
    counters[byteArr[7 + 8 * i]]++;
  }

  for (int i = 255; i > 127; i--) {
    counters[i] += acc;
    acc = counters[i];
  }

  for (int i = 0; i < 128; i++) {
    int tmp = counters[i];
    counters[i] = acc;
    acc += tmp;
  }

  for (int i = 0; i < size; i++) {
    auto byteIndex = offset + 8 * i;
    auto byte = byteArr[byteIndex];
    if (byte < 128) {
      result[counters[byte]] = source[i];
      counters[byte]++;
    } else {
      counters[byte]--;
      result[counters[byte]] = source[i];
    }
  }

  return result;
}

std::vector<double> serialRadixSort(const std::vector<double>& source) {
  auto tmp = std::vector<double>(source.size());
  auto tmpForSwap = std::vector<double>(source.size());
  auto copy = std::vector<double>(source);

  for (int i = 0; i < 8; i++) {
    tmp = CountingSort(copy, tmp, i);
    tmpForSwap = copy;
    copy = tmp;
    tmp = tmpForSwap;
  }

  return SignSort(copy);
}

std::vector<double> RadixSort(const std::vector<double>& source) {
  int threadCount = 4;
  std::vector<double> resultVector;
  std::vector<std::vector<double>> sortedVectors;

  int bunch_size = (source.size() - 1) / threadCount + 1;
  std::vector<std::vector<double>> bunches((source.size() + bunch_size) /
                                           bunch_size);
  int bunchesCount = bunches.size();

  for (size_t i = 0; i < source.size(); i += bunch_size) {
    auto last = std::min(source.size(), i + bunch_size);
    auto index = i / bunch_size;
    auto& vec = bunches[index];
    vec.reserve(last - i);
    move(source.begin() + i, source.begin() + last, back_inserter(vec));
  }

  bunchesCount = bunches.size();

#pragma omp parallel num_threads(threadCount)
  {
#pragma omp for
    for (auto i = 0; i < bunchesCount; i++) {
      auto sortedBunch = serialRadixSort(bunches[i]);
#pragma omp critical
      sortedVectors.push_back(sortedBunch);
      sortedBunch.clear();
    }
  }

  resultVector = sortedVectors[0];
  if (bunchesCount != 1) {
    for (auto i = 0; i < bunchesCount - 1; i++) {
      resultVector = Merge(resultVector, sortedVectors[i + 1]);
    }
  }
  sortedVectors.clear();
  bunches.clear();

  return resultVector;
}

std::vector<double> Merge(const std::vector<double>& left,
                          const std::vector<double>& right) {
  std::vector<double> result((left.size() + right.size()));

  auto const leftSize = static_cast<int>(left.size());
  auto const rightSize = static_cast<int>(right.size());

  auto i = 0, j = 0, k = 0;
  for (k = 0; k < (leftSize + rightSize - 1) && i < leftSize && j < rightSize;
       k++) {
    if (left[i] < right[j])
      result[k] = left[i++];
    else
      result[k] = right[j++];
  }

  while (i < leftSize) {
    result[k++] = left[i++];
  }
  while (j < rightSize) {
    result[k++] = right[j++];
  }

  return result;
}
