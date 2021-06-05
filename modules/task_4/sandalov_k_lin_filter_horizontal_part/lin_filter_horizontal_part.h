// Copyright 2020 Sandalov Konstantin
#ifndef MODULES_TASK_4_SANDALOV_K_LIN_FILTER_HORIZONTAL_PART_LIN_FILTER_HORIZONTAL_PART_H_
#define MODULES_TASK_4_SANDALOV_K_LIN_FILTER_HORIZONTAL_PART_LIN_FILTER_HORIZONTAL_PART_H_

#include <vector>
#include <cstdint>

namespace my {

template<typename type>
struct Image {
    type* image{nullptr};
    int64_t height{0};
    int64_t width{0};
    int8_t dims{0};
    Image(type* image, int64_t height, int64_t width, int8_t dims):
        image(image),
        height(height),
        width(width),
        dims(dims)
    {}
};


std::vector<float> createGaussianFilter(int64_t N, double sigma);
std::vector<uint8_t> linearFilterSeq(const Image<uint8_t>& srcImg,
    const std::vector<float>& kernel, int64_t K);
std::vector<uint8_t> linearFilterParStd(const Image<uint8_t>& srcImg,
    const std::vector<float>& kernel, int64_t K);
}  // namespace my

#endif  // MODULES_TASK_4_SANDALOV_K_LIN_FILTER_HORIZONTAL_PART_LIN_FILTER_HORIZONTAL_PART_H_
