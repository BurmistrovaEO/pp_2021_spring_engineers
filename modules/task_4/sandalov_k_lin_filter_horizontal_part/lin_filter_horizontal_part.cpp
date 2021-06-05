// Copyright 2020 Sandalov Konstantin
#include <unapproved.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include "../../modules/task_4/sandalov_k_lin_filter_horizontal_part/lin_filter_horizontal_part.h"

#define _USE_MATH_DEFINES

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

namespace my {

std::int64_t clamp(std::int64_t val, std::int64_t min, std::int64_t max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

std::vector<float> createGaussianFilter(std::int64_t N, double sigma) {
    if (N % 2 == 0) {
        throw std::runtime_error("Kernel size should be odd");
    }
    if (N <= 0) {
        throw std::runtime_error("Kernel size should be positive");
    }
    if (sigma <= 0) {
        throw std::runtime_error("Sigma should be positive");
    }

    std::vector<float> imgFilter(N*N);
    std::int64_t n = N / 2;
    float sum{};
    for (std::int64_t i = -n; i <= n; ++i) {
        for (std::int64_t j = -n; j <= n; ++j) {
            double r = (i * i + j * j);
            double s{ 2.0 * sigma * sigma };

            imgFilter[(n + i) * N + (n + j)] = static_cast<float>
                (exp(-(r / s)) / (M_PI * s));
            sum += imgFilter[(n + i) * N + (n + j)];
        }
    }
    for (std::int64_t i = 0; i < N; ++i) {
        for (std::int64_t j = 0; j < N; ++j) {
            imgFilter[i * N + j] /= sum;
        }
    }
    return imgFilter;
}

std::uint8_t applyKernel(const Image<std::uint8_t>& srcImg,
    const std::vector<float>& kernel, std::int64_t X,
                                    std::int64_t Y,
                                    std::int64_t K) {
    float resValue{0};
    std::int64_t width = K / 2;
    for (std::int64_t i = -width; i <= width; ++i) {
        for (std::int64_t j = -width; j <= width; ++j) {
            std::int64_t curX = clamp(X + i, 0, srcImg.height - 1);
            std::int64_t curY = clamp(Y + j * srcImg.dims, 0, srcImg.width*srcImg.dims - 1);
            resValue += (srcImg.image[(curX)*srcImg.width*srcImg.dims + (curY)]
                * kernel[(width + i) * K + (width + j)]);
        }
    }
    return static_cast<std::uint8_t>(resValue);
}

std::vector<std::uint8_t> linearFilterSeq(const Image<std::uint8_t>& srcImg,
    const std::vector<float>& kernel, std::int64_t K) {
    if (srcImg.height <= 0 || srcImg.width <= 0 || K <= 0) {
        throw std::runtime_error("Arguments should be positive");
    }
    if (kernel.size() != static_cast<std::size_t>(K*K)) {
        throw std::runtime_error("Incorrect Kernel Size");
    }
    std::vector<std::uint8_t> resImg(srcImg.height * srcImg.width * srcImg.dims);
    for (std::int64_t j = 0; j < srcImg.height; ++j) {
        for (std::int64_t f = 0; f < srcImg.width; ++f) {
            for (std::int8_t dim = 0; dim < srcImg.dims; ++dim)
                resImg[j * srcImg.width * srcImg.dims + f * srcImg.dims + dim] = applyKernel(srcImg,
                    kernel, j, f * srcImg.dims + dim, K);
        }
    }
    return resImg;
}

template<class Function>
bool my_parallel_for(const int64_t first, const int64_t last, Function&& f,
    const int num_threads = 1, const int64_t threshold = 1) {
    if (threshold < 1 || num_threads < 1) return false;
    const int64_t group_size = std::max(threshold, (last - first) / num_threads);
    std::vector<std::thread> my_threads;
    for (int64_t it = first; it < last; it += group_size) {
        if (it == last - 1) {
            my_threads.emplace_back(std::thread(f, it));
        } else {
            my_threads.emplace_back(std::thread([&f, it, last, group_size](){
                const auto cur_it = std::min(it, last);
                const auto cur_last = std::min(it + group_size, last);
                for (int64_t t = cur_it; t < cur_last; ++t) {
                    f(t);
                }
            }));
        }
    }
    std::for_each(my_threads.begin(), my_threads.end(), [](std::thread& thread) {
        thread.join();
    });
    return true;
}

std::vector<std::uint8_t> linearFilterParStd(const Image<std::uint8_t>& srcImg,
    const std::vector<float>& kernel, std::int64_t K) {
    if (srcImg.height <= 0 || srcImg.width <= 0 || K <= 0) {
        throw std::runtime_error("Arguments should be positive");
    }
    if (kernel.size() != static_cast<size_t>(K*K)) {
        throw std::runtime_error("Incorrect Kernel Size");
    }
    std::vector<std::uint8_t> resImg(srcImg.height * srcImg.width * srcImg.dims);

    my_parallel_for(int64_t{0}, srcImg.height, [&resImg, &srcImg, &kernel, K](int64_t index){
        for (std::int64_t f = 0; f < srcImg.width; ++f) {
            for (std::int8_t dim = 0; dim < srcImg.dims; ++dim)
                resImg[index * srcImg.width * srcImg.dims + f * srcImg.dims + dim] = applyKernel(srcImg,
                    kernel, index, f * srcImg.dims + dim, K);
        }
    }, 4);
    return resImg;
}
}  // namespace my
