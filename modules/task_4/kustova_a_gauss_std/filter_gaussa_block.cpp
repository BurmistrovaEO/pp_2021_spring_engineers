// Copyright 2021 Kustova Anastasiya
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include "../../../3rdparty/unapproved/unapproved.h"
#include "../../../modules/task_4/kustova_a_gauss_std/filter_gaussa_block.h"

std::vector<int> resultImage;

void threadFunction(std::vector<std::vector<int>> array, int i, int count_thread,
    int width, int height, const std::vector<int> & img, const std::vector<float> & kernel, int radius) {
    for (int t = i; t < static_cast<int>(array.size()); t++) {
        if ((t % count_thread) == i) {
            int j_start = array[t][0];
            int j_finish = array[t][0] + array[t][2];
            int i_start = array[t][1];
            int i_finish = array[t][1] + array[t][3];
            for (int i = i_start; i < i_finish && i < width; i++) {
                for (int j = j_start; j < j_finish && j < height; j++) {
                    int color = calculateNewPixelColor(img, width, height, i, j, radius, kernel);
                    resultImage[j * width + i] = color;
                }
            }
            t += count_thread - 1;
        }
    }
}

std::vector<int> gaussianFilter(const std::vector<int> & img, int width, int height,
    int radius, float sigma, int count_thread) {
    resultImage = img;
    int size = 2 * radius + 1;
    int l = 0;
    int k = 0;
    std::vector<float> kernel(size * size);
    kernel = createGaussianKernel(radius, sigma);
    std::vector<std::vector<int>> array;
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
            array.push_back(tup);
            k += block_w;
        }
        l += block_h;
    }
    std::vector<std::thread> threads;
    for (int i = 0; i < count_thread; i++) {
        threads.push_back(std::thread(threadFunction, array, i, count_thread,
            width, height, std::ref(img), std::ref(kernel), radius));
    }
    for (int i = 0; i < count_thread; i++) {
        threads[i].join();
    }
    return resultImage;
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
    double sumColor;
    int neighborColor;
    int idX, idY;
    int l = 0;
    int k = 0;
    sumColor = 0;
    for (l = -radius; l <= radius; l++) {
        for (k = -radius; k <= radius; k++) {
            idX = Clamp(x + k, 0, width - 1);
            idY = Clamp(y + l, 0, height - 1);
            neighborColor = img[idY * width + idX];
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

