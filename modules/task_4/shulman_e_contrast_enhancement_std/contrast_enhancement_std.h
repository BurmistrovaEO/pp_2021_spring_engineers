// Copyright 2021 Shulman Egor
#ifndef MODULES_TASK_4_SHULMAN_E_CONTRAST_ENHANCEMENT_STD_CONTRAST_ENHANCEMENT_STD_H_
#define MODULES_TASK_4_SHULMAN_E_CONTRAST_ENHANCEMENT_STD_CONTRAST_ENHANCEMENT_STD_H_

#include <algorithm>
#include <vector>

std::vector<int> getRandomMatrix(int n, int m);
std::vector<int> ContrastEnhancementSeq(const std::vector<int>& matrix);
std::vector<int> ContrastEnhancementStd(const std::vector<int>& matrix);

#endif  // MODULES_TASK_4_SHULMAN_E_CONTRAST_ENHANCEMENT_STD_CONTRAST_ENHANCEMENT_STD_H_
