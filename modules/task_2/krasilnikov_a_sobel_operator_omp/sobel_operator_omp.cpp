// Copyright 2021 Krasilnikov Alexey
#include "../../modules/task_2/krasilnikov_a_sobel_operator_omp/sobel_operator_omp.h"

#include <cmath>
#include <algorithm>
#include <vector>

void SobelOperator_parallel(unsigned char* image, int height, int weight) {
    std::vector<unsigned char> ans(height * weight);
#pragma omp parallel for
    for (int i = 1; i < height - 1; ++i) {
        for (int j = 1; j < weight - 1; ++j) {
            auto Gy = image[(i + 1) * weight + (j - 1)] -
                      image[(i - 1) * weight + (j - 1)] +
                      2 * image[(i + 1) * weight + (j)] -
                      2 * image[(i - 1) * weight + (j)] +
                      image[(i + 1) * weight + (j + 1)] -
                      image[(i - 1) * weight + (j + 1)];
            auto Gx = image[(i + 1) * weight + (j + 1)] -
                      image[(i - 1) * weight + (j - 1)] +
                      2 * image[(i) * weight + (j + 1)] -
                      2 * image[(i) * weight + (j - 1)] +
                      image[(i - 1) * weight + (j + 1)] -
                      image[(i + 1) * weight + (j - 1)];
            auto pixel = round(
                    std::max(0., std::min(255., sqrt(Gy * Gy + Gx * Gx))));
            ans[i * weight + j] = static_cast<unsigned char>(pixel);
        }
    }
#pragma omp parallel for
    for (int i = 0; i < height * weight; ++i) {
        image[i] = ans[i];
    }
}

void SobelOperator(unsigned char* image, int height, int weight) {
    std::vector<unsigned char> ans(height * weight);
    for (int i = 1; i < height - 1; ++i) {
        for (int j = 1; j < weight - 1; ++j) {
            auto Gy = image[(i + 1) * weight + (j - 1)] -
                      image[(i - 1) * weight + (j - 1)] +
                      2 * image[(i + 1) * weight + (j)] -
                      2 * image[(i - 1) * weight + (j)] +
                      image[(i + 1) * weight + (j + 1)] -
                      image[(i - 1) * weight + (j + 1)];
            auto Gx = image[(i + 1) * weight + (j + 1)] -
                      image[(i - 1) * weight + (j - 1)] +
                      2 * image[(i) * weight + (j + 1)] -
                      2 * image[(i) * weight + (j - 1)] +
                      image[(i - 1) * weight + (j + 1)] -
                      image[(i + 1) * weight + (j - 1)];
            auto pixel = round(
                    std::max(0., std::min(255., sqrt(Gy * Gy + Gx * Gx))));
            ans[i * weight + j] = static_cast<unsigned char>(pixel);
        }
    }
    for (int i = 0; i < height * weight; ++i) {
        image[i] = ans[i];
    }
}
