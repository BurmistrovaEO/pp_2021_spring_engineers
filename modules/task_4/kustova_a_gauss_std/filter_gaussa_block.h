// Copyright 2021 Kustova Anastasiya
#ifndef MODULES_TASK_4_KUSTOVA_A_GAUSS_STD_FILTER_GAUSSA_BLOCK_H_
#define MODULES_TASK_4_KUSTOVA_A_GAUSS_STD_FILTER_GAUSSA_BLOCK_H_
#include <vector>

std::vector<int> gaussianFilter(const std::vector<int> & img, int width,
    int height, int radius, float sigma, int count_thread);
std::vector<float> createGaussianKernel(int radius, float sigma);
int calculateNewPixelColor(std::vector<int> img, int width, int height,
    int x, int y, int radius, std::vector<float> kernel);
int Clamp(int value, int min, int max);
void threadFunction(std::vector<std::vector<int>> array, int i, int count_thread,
    int width, int height, const std::vector<int> & img, const std::vector<float> & kernel, int radius);
#endif  // MODULES_TASK_4_KUSTOVA_A_GAUSS_STD_FILTER_GAUSSA_BLOCK_H_
