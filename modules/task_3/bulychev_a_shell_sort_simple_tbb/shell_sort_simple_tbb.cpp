// Copyright 2021 Bulychev Andrey
#include <tbb/tbb.h>
#include <omp.h>
#include <random>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>
#include "../../../modules/task_3/bulychev_a_shell_sort_simple_tbb/shell_sort_simple_tbb.h"


std::vector<double> getRandomVector(int size) {
  if (size < 0) {
    throw "Error size";
  }
  std::vector<double> vec(size);
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  for (int i = 0; i < size; i++) vec[i] = gen() % 10000;
  return vec;
}

std::vector<double> Shell_sort(const std::vector<double>& vec) {
  std::vector<double> tmp(vec);
  for (int step = tmp.size() / 2; step != 0; step /= 2) {
    for (int i = step; i < static_cast<int>(tmp.size()); i++) {
      for (int j = i - step; j >= 0 && tmp[j] > tmp[j + step]; j -= step) {
        double aV = tmp[j];
        tmp[j] = tmp[j + step];
        tmp[j + step] = aV;
      }
    }
  }
  return tmp;
}

std::vector<double> Merge(const std::vector<double>& vec1,
                          const std::vector<double>& vec2) {
  int size1 = vec1.size();
  int size2 = vec2.size();
  std::vector<double> tmp(size1 + size2);
  int i = 0, j = 0;
  for (int k = 0; k < size1 + size2; k++) {
    if (i > size1 - 1) {
      double a = vec2[j];
      tmp[k] = a;
      j++;
    } else {
      if (j > size2 - 1) {
        double a = vec1[i];
        tmp[k] = a;
        i++;
      } else {
        if (vec1[i] < vec2[j]) {
          double a = vec1[i];
          tmp[k] = a;
          i++;
        } else {
          double b = vec2[j];
          tmp[k] = b;
          j++;
        }
      }
    }
  }
  return tmp;
}


std::vector<double> Shell_Sort_tbb(const std::vector<double>& vec, int num_threads) {
    std::vector<double> copy(vec);
    if (vec.size() == 1) {
        return copy;
    }

    if (static_cast<int>(vec.size()) < num_threads) {
      return Shell_sort(vec);
    }
    int local_size = vec.size() / num_threads;

    std::vector<double> tmp;
    std::vector<std::vector<double>> local_vec;

    for (int num_thread = 0; num_thread < num_threads; num_thread++) {
        if (num_thread == num_threads - 1) {
            tmp.insert(tmp.end(), copy.begin() + local_size * num_thread, copy.end());
        } else {
            tmp.insert(tmp.end(), copy.begin() + local_size * num_thread,
                copy.begin() + local_size * (num_thread + 1));
        }
        local_vec.push_back(tmp);
        tmp.clear();
    }

    tbb::task_scheduler_init init(num_threads);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, local_vec.size(), 1),
    [&local_vec](const tbb::blocked_range<size_t>& range){
        for (size_t i = range.begin(); i != range.end(); ++i) {
            local_vec[i] = Shell_sort(local_vec[i]);
        }
    }, tbb::simple_partitioner());

    init.terminate();
    std::vector<double> merged = local_vec[0];
    for (size_t i = 1; i < local_vec.size(); ++i) {
        merged = Merge(merged, local_vec[i]);
    }

    return merged;
}
