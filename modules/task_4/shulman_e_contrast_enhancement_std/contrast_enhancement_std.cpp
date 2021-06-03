// Copyright 2021 Shulman Egor
#include <time.h>
#include <random>
#include <iostream>
#include "../../../modules/task_4/shulman_e_contrast_enhancement_std/contrast_enhancement_std.h"
#include "../../../3rdparty/unapproved/unapproved.h"

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

std::vector<int> ContrastEnhancementStd(const std::vector<int>& matrix) {
    std::vector<int> result(matrix);

    auto temp = std::minmax_element(matrix.begin(), matrix.end());
    int yMax = *(temp.second);
    int yMin = *(temp.first);

    auto threadFill = [&](std::vector<int>::iterator first,
                          std::vector<int>::iterator last) {
        for (auto it = first; it != last; ++it) {
            *it = ((*it - yMin) * 255) / (yMax - yMin);
        }
    };

    const size_t size = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(size);
    auto block = result.size() / size;
    auto work_iter = std::begin(result);

    for (auto it = std::begin(threads); it != std::end(threads) - 1; ++it) {
        *it = std::thread(threadFill, work_iter, work_iter + block);
        work_iter += block;
    }
    threads.back() = std::thread(threadFill, work_iter, std::end(result));

    for (auto&& i : threads) {
        i.join();
    }
    return result;
}
