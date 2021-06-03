// Copyright 2021 Ivanov Yaroslav
#include <omp.h>
#include <random>
#include <ctime>
#include <algorithm>
#include <vector>
#include <iostream>
#include "../../../modules/task_2/ivanov_y_double_sort_opm/sort_double.h"

std::vector<double> BubbleForCompare(std::vector<double> inputNumbers, const int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (inputNumbers[i] < inputNumbers[j]) {
                double tmp = inputNumbers[i];
                inputNumbers[i] = inputNumbers[j];
                inputNumbers[j] = tmp;
            }
        }
    }
    return inputNumbers;
}

std::vector<double> generateDouble(const int n) {
    std::mt19937 gen;
    gen.seed(static_cast<int>(time(0)));
    std::vector<double> vec(n);
    for (int i = 0; i < n; i++) {
        double numb = gen() / 1000;
        while (numb >= 1001 || numb <= -1001) {
            numb /= 100;
        }
        vec[i] = numb;
    }
    return vec;
}

std::vector<double> transit_omp(std::vector<double> inputNumbers, std::vector<double> loc, const int size, int add) {
    unsigned char* pmem = reinterpret_cast<unsigned char*>(inputNumbers.data());
    int counters[256] = { 0 };
    int sum = 0;
    // omp_set_num_threads(1);
    #pragma omp for
    for (int i = 0; i < size; i++) {
        counters[pmem[8 * i + add]]++;
    }
    #pragma omp for
    for (int i = 0; i < 256; i++) {
        int tmp = counters[i];
        counters[i] = sum;
        sum += tmp;
    }
    #pragma omp for
    for (int i = 0; i < size; i++) {
        int index = 8 * i + add;
        loc[counters[pmem[index]]] = inputNumbers[i];
        counters[pmem[index]]++;
    }
    return loc;
}

std::vector<double> loc_sort_omp(std::vector<double> inputNumbers, const int size) {
    std::vector<double> loc(size);
    std::vector<double> tmp(size);
    #pragma omp for
    for (int i = 0; i < size; i++) {
        loc[i] = 0;
        tmp[i] = 0;
    }
    for (int i = 0; i < 8; i++) {
        loc = transit_omp(inputNumbers, loc, size, i);
        tmp = inputNumbers;
        inputNumbers = loc;
        loc = tmp;
    }
    return inputNumbers;
}

std::vector<double> transit(std::vector<double> inputNumbers, std::vector<double> loc, const int size, int add) {
    unsigned char* pmem = reinterpret_cast<unsigned char*>(inputNumbers.data());
    int counters[256] = { 0 };
    int sum = 0;

    for (int i = 0; i < size; i++) {
        counters[pmem[8 * i + add]]++;
    }
    for (int i = 0; i < 256; i++) {
        int tmp = counters[i];
        counters[i] = sum;
        sum += tmp;
    }
    for (int i = 0; i < size; i++) {
        int index = 8 * i + add;
        loc[counters[pmem[index]]] = inputNumbers[i];
        counters[pmem[index]]++;
    }
    return loc;
}

std::vector<double> loc_sort(std::vector<double> inputNumbers, const int size) {
    std::vector<double> loc(size);
    std::vector<double> tmp(size);
    for (int i = 0; i < size; i++) {
        loc[i] = 0;
        tmp[i] = 0;
    }
    for (int i = 0; i < 8; i++) {
        loc = transit(inputNumbers, loc, size, i);
        tmp = inputNumbers;
        inputNumbers = loc;
        loc = tmp;
    }
    return inputNumbers;
}
