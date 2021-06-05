// Copyright 2021 Solomakhin Sergey
#include "../../../modules/task_3/solomakhin_s_hoar_sort_simple_merge/hoar_sort_simple_merge.h"

#include <tbb/tbb.h>
#include <random>
#include <ctime>
#include <vector>
#include <climits>
#include <functional>
#include <algorithm>
#include <stdexcept>

std::vector<int> random_gen(int size, int type) {
    std::vector<int> newVec(size);

    std::random_device dev;
    std::mt19937 randomGenerator(dev());

    if (type == 1) {
        std::uniform_int_distribution<int> randomValue(1, INT_MAX);

        for (int& item : newVec) {
            item = randomValue(randomGenerator);
        }
    } else if (type >= 0 && type <= 3) {
        std::uniform_int_distribution<int> randomValue(INT_MIN, INT_MAX);

        for (int& item : newVec) {
            item = randomValue(randomGenerator);
        }

        if (type == 2) {
            sort(newVec.begin(), newVec.end());
        } else if (type == 3) {
            sort(newVec.begin(), newVec.end(), std::greater<int>());
        }
    } else {
        throw std::invalid_argument("Invalid type of vector chosen");
    }

    return newVec;
}

int partition(std::vector<int>* arr, int left, int right) {
    int pivot = arr->at(right);
    int smallerElementIdx = left - 1;

    for (int i = left; i < right; ++i) {
        if (arr->at(i) < pivot) {
            std::swap(arr->at(++smallerElementIdx), arr->at(i));
        }
    }

    std::swap(arr->at(smallerElementIdx + 1), arr->at(right));

    return smallerElementIdx + 1;
}

void hoar_sort(std::vector<int>* arr, int left, int right) {
    if (left < right) {
        int partitionIdx = partition(arr, left, right);

        hoar_sort(arr, left, partitionIdx - 1);
        hoar_sort(arr, partitionIdx + 1, right);
    }
}

std::vector<int>merge(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> result((a.size() + b.size()));

    int i = 0, j = 0, k = 0;
    while (i < static_cast<int>(a.size()) && j < static_cast<int>(b.size())) {
        if (a[i] < b[j])
            result[k] = a[i++];
        else
            result[k] = b[j++];
        k++;
    }
    while (i < static_cast<int>(a.size())) {
        result[k++] = a[i++];
    }
    while (j < static_cast<int>(b.size())) {
        result[k++] = b[j++];
    }

    return result;
}

void hoar_sort_tbb(std::vector<int>* arr, int numThreads) {
    if (static_cast<int>(arr->size()) < numThreads) {
        hoar_sort(arr, 0, arr->size() - 1);
        return;
    }
    if (arr->size() == 1) {
        return;
    }
    tbb::task_scheduler_init init(numThreads);

    int elementsPerSegment = arr->size() / numThreads;

    std::vector<std::vector<int>> segments;
    std::vector<int> segment;

    for (int i = 0; i < numThreads; ++i) {
        if (i == numThreads - 1) {
            segment.insert(segment.end(),
                arr->begin() + i * elementsPerSegment,
                arr->end());
        } else {
            segment.insert(segment.end(),
                arr->begin() + i * elementsPerSegment,
                arr->begin() + (i + 1) * elementsPerSegment);
        }

        segments.push_back(segment);
        segment.clear();
    }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, segments.size(), 1),
        [&segments](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                hoar_sort(&segments[i], 0, segments[i].size() - 1);
            }
        }, tbb::simple_partitioner());

    *arr = segments[0];
    for (size_t i = 1; i < segments.size(); ++i) {
        *arr = merge(*arr, segments[i]);
    }
}
