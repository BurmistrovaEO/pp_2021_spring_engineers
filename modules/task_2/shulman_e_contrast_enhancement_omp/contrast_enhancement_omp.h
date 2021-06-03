// Copyright 2021 Shulman Egor
#ifndef MODULES_TASK_2_SHULMAN_E_CONTRAST_ENHANCEMENT_OMP_CONTRAST_ENHANCEMENT_OMP_H_
#define MODULES_TASK_2_SHULMAN_E_CONTRAST_ENHANCEMENT_OMP_CONTRAST_ENHANCEMENT_OMP_H_

#include <vector>

std::vector<int> getRandomMatrix(int n, int m);
std::vector<int> ContrastEnhancement(const std::vector<int>& matrix);
std::vector<int> ContrastEnhancementOMP(const std::vector<int>& matrix);

#endif  // MODULES_TASK_2_SHULMAN_E_CONTRAST_ENHANCEMENT_OMP_CONTRAST_ENHANCEMENT_OMP_H_
