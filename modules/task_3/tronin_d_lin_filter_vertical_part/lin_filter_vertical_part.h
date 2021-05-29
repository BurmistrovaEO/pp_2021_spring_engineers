// Copyright 2021 Tronin Dmitry
#ifndef MODULES_TASK_3_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_
#define MODULES_TASK_3_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_

#include <tbb/tbb.h>
#include <cstdint>
#include <vector>

std::vector<double> CalculateGaussFilterTBB(size_t size, double sigma);
std::vector<uint8_t> ApplyFilterTBB(const std::vector<double> &filter,
                                    const std::vector<uint8_t> &image,
                                    size_t width,
                                    size_t height,
                                    size_t filter_size, size_t number_of_threads);

template<class T>
inline T ClampTBB(T value, T low, T high) {
    if (value < low) {
        return low;
    }
    if (value > high) {
        return high;
    }
    return value;
}

#endif  // MODULES_TASK_3_TRONIN_D_LIN_FILTER_VERTICAL_PART_LIN_FILTER_VERTICAL_PART_H_
