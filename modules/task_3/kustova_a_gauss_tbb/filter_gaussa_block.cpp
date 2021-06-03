// Copyright 2021 Kustova Anastasiya
#include <math.h>
#include <tbb/tbb.h>
#include <iostream>
#include <vector>
#include <tuple>
#include "../../../modules/task_3/kustova_a_gauss_tbb/filter_gaussa_block.h"

class GaussianParallel {
 private:
    const std::vector<int> &img;
    std::vector<int> *res;
    const std::vector<float> &kernel;
    const std::vector<std::vector<int>> &array;
    int width, height, radius, block_width, block_height;

 public:
    void operator()(const tbb::blocked_range<int> &r) const {
        for (int i = r.begin(); i != r.end(); ++i) {
            int j_start = array[i][0];
            int j_finish = array[i][0] + block_height;
            int i_start = array[i][1];
            int i_finish = array[i][1] + block_width;
            for (int i = i_start; i < i_finish && i < width; i++) {
                for (int j = j_start; j < j_finish && j < height; j++) {
                    int color = calculateNewPixelColor(img, width, height, i, j, radius, kernel);
                    res->at(j * width + i) = color;
                }
            }
        }
    }

    GaussianParallel(const std::vector<int> &img, std::vector<int> *res, const std::vector<float> &kernel,
        const std::vector<std::vector<int>> &array, int width, int height, int radius,
        int block_width, int block_height) :
        img(img),
        res(res),
        kernel(kernel),
        array(array),
        width(width),
        height(height),
        radius(radius),
        block_width(block_width),
        block_height(block_height){}
};

std::vector<int> gaussianFilter(const std::vector<int> & img, int width,
    int height, int radius, float sigma, int count_thread) {
    std::vector<int> res(img);
    int size = 2 * radius + 1;
    std::vector<float> kernel(size * size);
    int l = 0;
    int k = 0;
    kernel = createGaussianKernel(radius, sigma);
    std::vector<std::vector<int>> array1;
    int grid_size = ceil(static_cast<double>(sqrt(count_thread)));
    int block_height = height / grid_size;
    int block_width = width / grid_size;
    while (l < height) {
        k = 0;
        int block_w = block_width;
        int block_h = block_height;
        while (k < width) {
            if (width - k < block_width + block_width / 2) {
                block_w = width - k;
            } else {
                block_w = block_width;
            }
            if (height - l < block_height / 2) {
                block_h = height - l;
            } else {
                block_h = block_height;
            }
            std::vector<int> tup = { l, k , block_h, block_w };
            array1.push_back(tup);
            k += block_w;
        }
        l += block_h;
    }
    tbb::task_scheduler_init init(count_thread);
    GaussianParallel tmp(img, &res, kernel, array1, width, height, radius, block_width, block_height);
    tbb::parallel_for(tbb::blocked_range<int>(0, static_cast<int>(array1.size())), tmp);
    return res;
}
std::vector<float> createGaussianKernel(int radius, float sigma) {
    int size = 2 * radius + 1;
    std::vector<float> kernel(size * size);
    float norm = 0;
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            kernel[(i + radius)* size + j + radius] = static_cast<float>(exp(-(i * i + j * j) / (2 * sigma * sigma)));
            norm = norm + kernel[(i + radius) * size + j + radius];
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            kernel[i * size + j] = kernel[i * size + j] / norm;
        }
    }
    return kernel;
}

int calculateNewPixelColor(std::vector<int> img, int width, int height,
    int x, int y, int radius, std::vector<float> kernel) {
    int size = radius * 2 + 1;
    double sumColor = 0;
    for (int l = -radius; l <= radius; l++) {
        for (int k = -radius; k <= radius; k++) {
            int idX = Clamp(x + k, 0, width - 1);
            int idY = Clamp(y + l, 0, height - 1);
            int neighborColor = img[idY * width + idX];
            sumColor += neighborColor * kernel[(k + radius) * size + l + radius];
        }
    }
    return Clamp(static_cast<int>(sumColor), 0, 255);
}
int Clamp(int value, int min, int max) {
    if (value < min) {
        return min;
    }
    if (value > max) {
        return max;
    }
    return value;
}

