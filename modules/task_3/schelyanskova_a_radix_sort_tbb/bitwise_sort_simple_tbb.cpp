// Copyright 2021 Schelyanskova Anastasiia
#include <tbb/tbb.h>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include "../../../modules/task_3/schelyanskova_a_radix_sort_tbb/bitwise_sort_simple_tbb.h"


std::vector<int> getRandomVector(int sz) {
    std::random_device dev;
    std::mt19937 gen(dev());
    std::vector<int> vect(sz);
    for (int i = 0; i < sz; i++) {
        vect[i] = gen() % 100;
    }
    return vect;
}

bool vector_sort(std::vector<int> vect) {
    int size = vect.size();
    for (int i = 0; i < size - 1; i++) {
        if (vect[i] > vect[i + 1]) {
            return false;
        }
    }
    return true;
}

int get_max_power(std::vector<int> vect) {
    int max_power = 0;
    int tmp = vect[0];
    int size = vect.size();
    for (int i = 0; i < size; i++) {
        if (tmp < vect[i]) {
            tmp = vect[i];
        }
    }
    if (tmp == 0) {
        return 1;
    }
    while (tmp != 0) {
        tmp = tmp / 10;
        max_power++;
    }
    return max_power;
}

std::vector<int> RadixSort(std::vector<int> vect, int size) {
    std::vector<std::vector<int>> vect_start(10);
    std::vector<int> vect1 = vect;
    int tmp = 0;
    int power = 1;
    int max_power;
    if (vector_sort(vect)) {
        return vect;
    }
    max_power = get_max_power(vect);
    while (power <= max_power) {
        for (int i = 0; i < size; i++) {
            tmp = (vect1[i] % static_cast<unsigned int>(pow(10, power))
                    / static_cast<unsigned int>(pow(10, power - 1)));
            vect_start[tmp].push_back(vect1[i]);
        }
        vect1.clear();
        for (int i = 0; i < 10; i++) {
            if (!vect_start[i].empty()) {
                int size_tmp = vect_start[i].size();
                for (int j = 0; j < size_tmp; j++) {
                    vect1.push_back(vect_start[i][j]);
                }
                vect_start[i].clear();
            }
        }
        power++;
    }
    return vect1;
}

std::vector<int> RadixSort_tbb(const std::vector<int>& v) {
    int num_threads = 4;
    if (static_cast<int>(v.size()) < num_threads) {
        return RadixSort(v, v.size());
    }
    std::vector<int> copy(v);
    if (v.size() == 1) {
        return copy;
    }

    int delta = v.size() / num_threads;

    std::vector<std::vector<int>> split_vec(0);
    std::vector<int> tmp;

    for (int thread_num = 0; thread_num < num_threads; thread_num++) {
        if (thread_num == num_threads - 1) {
            tmp.insert(tmp.end(), copy.begin() + delta * thread_num, copy.end());
        } else {
            tmp.insert(tmp.end(), copy.begin() + delta * thread_num,
                       copy.begin() + delta * (thread_num + 1));
        }
        split_vec.push_back(tmp);
        tmp.clear();
    }

    tbb::task_scheduler_init init(num_threads);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, split_vec.size(), 1),
                      [&split_vec](const tbb::blocked_range<size_t>& range){
                          for (size_t i = range.begin(); i != range.end(); ++i) {
                              split_vec[i] = RadixSort(split_vec[i], split_vec[i].size());
                          }
                      }, tbb::simple_partitioner());

    std::vector<int> merged = split_vec[0];
    for (size_t i = 1; i < split_vec.size(); ++i) {
        merged = Merge(merged, split_vec[i]);
    }

    return merged;
}

std::vector<int> Merge(const std::vector<int>& vect_left,
                       const std::vector<int>& vect_right) {
    std::vector<int> result((vect_left.size() + vect_right.size()));

    int i = 0, j = 0, k = 0;
    while (i < static_cast<int>(vect_left.size()) && j < static_cast<int>(vect_right.size())) {
        if (vect_left[i] < vect_right[j])
            result[k] = vect_left[i++];
        else
            result[k] = vect_right[j++];
        k++;
    }
    while (i < static_cast<int>(vect_left.size())) {
        result[k++] = vect_left[i++];
    }
    while (j < static_cast<int>(vect_right.size())) {
        result[k++] = vect_right[j++];
    }

    return result;
}
