// Copyright 2021 Krasilnikov Alexey
#include "../../modules/task_3/krasilnikov_a_sobel_operator_tbb/sobel_operator_tbb.h"

#include <tbb/tbb.h>
#include <cmath>
#include <vector>

void SobelOperator_parallel(unsigned char* image, int height, int weight) {
    std::vector<unsigned char> ans(height * weight);
    tbb::parallel_for(tbb::blocked_range<int> (1, height - 1, 4),
                      [&image, &ans, weight](
                          tbb::blocked_range<int> r) {
    for (int i = r.begin(); i < r.end(); ++i) {
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
            auto pp = sqrt(Gy * Gy + Gx * Gx);
            if (pp < 0) {
                pp = 0;
            }
            if (pp > 255) {
                pp = 255;
            }
            auto pixel = round(pp);
            ans[i * weight + j] = static_cast<unsigned char>(pixel);
        }
    }});
    tbb::parallel_for(tbb::blocked_range<int> (0, height * weight, 4),
                      [&image, &ans](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i) {
            image[i] = ans[i];
        }
    });
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
            auto pp = sqrt(Gy * Gy + Gx * Gx);
            if (pp < 0) {
                pp = 0;
            }
            if (pp > 255) {
                pp = 255;
            }
            auto pixel = round(pp);
            ans[i * weight + j] = static_cast<unsigned char>(pixel);
        }
    }
    for (int i = 0; i < height * weight; ++i) {
        image[i] = ans[i];
    }
}
