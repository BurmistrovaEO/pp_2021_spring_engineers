// Copyright 2021 Vizgalov Anton

#include <tbb/tbb.h>
#include <random>
#include <limits>
#include <functional>
#include <algorithm>
#include <stdexcept>

#include "../../../modules/task_3/vizgalov_a_hoare_sort_simple_merge/hoare_sort_simple_merge.h"

int partition(std::vector<int>* vec, int low, int high) {
    int pivot = vec->at(high);
    int smallerElementIdx = low - 1;

    for (int i = low; i < high; ++i) {
        if (vec->at(i) < pivot) {
            std::swap(vec->at(++smallerElementIdx), vec->at(i));
        }
    }

    std::swap(vec->at(smallerElementIdx + 1), vec->at(high));

    return smallerElementIdx + 1;
}

void hoareSort(std::vector<int>* vec, int low, int high) {
    if (low < high) {
        int partitionIdx = partition(vec, low, high);

        hoareSort(vec, low, partitionIdx - 1);
        hoareSort(vec, partitionIdx + 1, high);
    }
}

std::vector<int> merge(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    std::vector<int> merged(vec1.size() + vec2.size());

    int i = 0;
    int j = 0;
    int k = 0;
    while (i < static_cast<int>(vec1.size()) && j < static_cast<int>(vec2.size())) {
        if (vec1[i] < vec2[j]) {
            merged[k] = vec1[i];
            i++;
        } else {
            merged[k] = vec2[j];
            j++;
        }
        k++;
    }
    while (i < static_cast<int>(vec1.size())) {
        merged[k] = vec1[i];
        k++;
        i++;
    }
    while (j < static_cast<int>(vec2.size())) {
        merged[k] = vec2[j];
        k++;
        j++;
    }

    return merged;
}

void runHoareSort(std::vector<int>* vec) {
    int numThreads = 1;
    if (static_cast<int>(vec->size()) < numThreads) {
        hoareSort(vec, 0, vec->size() - 1);
        return;
    }
    if (vec->size() == 1) {
        return;
    }
    tbb::task_scheduler_init init(numThreads);

    int elementsPerSegment = vec->size() / numThreads;

    std::vector<std::vector<int>> segments;
    std::vector<int> segment;

    for (int i = 0; i < numThreads; ++i) {
        if (i == numThreads - 1) {
            segment.insert(segment.end(), vec->begin() + elementsPerSegment * i, vec->end());
        } else {
            segment.insert(segment.end(), vec->begin() + elementsPerSegment * i,
                vec->begin() + elementsPerSegment * (i + 1));
        }
        segments.push_back(segment);
        segment.clear();
    }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, segments.size(), 1),
    [&segments](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
            hoareSort(&segments[i], 0, segments[i].size() - 1);
        }
    }, tbb::simple_partitioner());

    std::vector<int> merged = segments[0];
    for (size_t i = 1; i < segments.size(); ++i) {
        merged = merge(merged, segments[i]);
    }

    *vec = merged;
}

// Types:
// 0 - any integers
// 1 - only positive integers
// 2 - any integers in non-decreasing order (sorted)
// 3 - any integers in non-increasing order (reverse-sorted)
std::vector<int> getRandomVector(int size, int type) {
    std::vector<int> newVec(size);

    std::random_device dev;
    std::mt19937 randomGenerator(dev());

    if (type == 1) {
        std::uniform_int_distribution<int> randomValue(1, std::numeric_limits<int>::max());

        for (int& item : newVec) {
            item = randomValue(randomGenerator);
        }
    } else if (type >= 0 && type <= 3) {
        std::uniform_int_distribution<int> randomValue(std::numeric_limits<int>::min(),
                                                       std::numeric_limits<int>::max());

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
