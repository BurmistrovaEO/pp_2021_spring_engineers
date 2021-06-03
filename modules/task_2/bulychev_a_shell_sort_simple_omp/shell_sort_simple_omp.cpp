// Copyright 2021 Bulychev Andrey
#include <omp.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>
#include "../../../modules/task_2/bulychev_a_shell_sort_simple_omp/shell_sort_simple_omp.h"


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


std::vector<double> Shell_sort_OMP(const std::vector<double>& vec, int num_threads) {
    std::vector<double> myvec(vec);
    if (vec.size() == 1) {
      return myvec;
    }
    if (static_cast<int>(vec.size()) < num_threads) {
        return Shell_sort(vec);
    }
    omp_set_num_threads(num_threads);

    int delta = vec.size() / num_threads;
    std::vector<double> local_vec;

#pragma omp parallel private(local_vec) shared(myvec, delta)
    {
        int num_thread = omp_get_thread_num();
        if (num_thread == num_threads - 1) {
            local_vec = std::vector<double>(myvec.size() - delta * (num_threads - 1));
        } else {
            local_vec = std::vector<double>(delta);
        }

        if (num_thread == num_threads - 1) {
            std::copy(myvec.begin() + delta * (num_threads - 1), myvec.end(), local_vec.begin());
        } else {
            std::copy(myvec.begin() + delta * num_thread,
                myvec.begin() + delta * (num_thread + 1), local_vec.begin());
        }
        local_vec = Shell_sort(local_vec);

#pragma omp barrier
        if (num_thread == 0) {
            myvec = local_vec;
        }

#pragma omp barrier
#pragma omp critical
        if (num_thread != 0) {
            myvec = Merge(myvec, local_vec);
        }
    }
    return myvec;
}
