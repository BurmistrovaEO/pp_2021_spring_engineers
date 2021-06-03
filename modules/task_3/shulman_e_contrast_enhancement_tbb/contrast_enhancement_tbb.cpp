// Copyright 2021 Shulman Egor
#include "../../../modules/task_3/shulman_e_contrast_enhancement_tbb/contrast_enhancement_tbb.h"
#include <time.h>
#include <random>
#include "tbb/parallel_for.h"

std::vector<int> getRandomMatrix(int n, int m) {
    if (n < 0 || m < 0)
        throw "Rows or columns less than 0";
    std::mt19937 random;
    random.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> matrix(n * m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i * n + j] = random() % 256;
        }
    }
    return matrix;
}

std::vector<int> ContrastEnhancementSeq(const std::vector<int>& matrix) {
    std::vector<int> result;
    int yMax = 0;
    int yMin = 255;
    for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
        yMax = yMax < matrix[i] ? matrix[i] : yMax;
        yMin = yMin > matrix[i] ? matrix[i] : yMin;
    }
    for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
        result.push_back(((matrix[i] - yMin) * 255) / (yMax - yMin));
    }
    return result;
}

std::vector<int> ContrastEnhancementOmp(const std::vector<int>& matrix) {
    std::vector<int> result(matrix.size());
    int yMax = 0;
    int yMin = 255;

#pragma omp parallel
    {
        int localMax = 0;
        int localMin = 255;
#pragma omp for
        for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
            localMax = localMax < matrix[i] ? matrix[i] : localMax;
            localMin = localMin > matrix[i] ? matrix[i] : localMin;
        }
#pragma omp critical
        {
            yMax = yMax > localMax ? yMax : localMax;
            yMin = yMin < localMin ? yMin : localMin;
        }
#pragma omp barrier
#pragma omp for
        for (int i = 0; i < static_cast<int>(matrix.size()); i++) {
            result[i] = ((matrix[i] - yMin) * 255) / (yMax - yMin);
        }
    }
    return result;
}

std::vector<int> ContrastEnhancementTbb(const std::vector<int>& matrix) {
    std::vector<int> result(matrix.size());
    auto temp = std::minmax_element(matrix.begin(), matrix.end());
    int yMax = *(temp.second);
    int yMin = *(temp.first);
    tbb::parallel_for(0, static_cast<int>(matrix.size()), [&](const int i) {
        result[i] = ((matrix[i] - yMin) * 255) / (yMax - yMin);
    });
    return result;
}
