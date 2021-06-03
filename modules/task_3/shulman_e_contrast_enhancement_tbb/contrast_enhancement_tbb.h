// Copyright 2021 Shulman Egor
#ifndef MODULES_TASK_3_SHULMAN_E_CONTRAST_ENHANCEMENT_TBB_CONTRAST_ENHANCEMENT_TBB_H_
#define MODULES_TASK_3_SHULMAN_E_CONTRAST_ENHANCEMENT_TBB_CONTRAST_ENHANCEMENT_TBB_H_

#include <algorithm>
#include <vector>

std::vector<int> getRandomMatrix(int n, int m);
std::vector<int> ContrastEnhancementSeq(const std::vector<int>& matrix);
std::vector<int> ContrastEnhancementOmp(const std::vector<int>& matrix);
std::vector<int> ContrastEnhancementTbb(const std::vector<int>& matrix);

#endif  // MODULES_TASK_3_SHULMAN_E_CONTRAST_ENHANCEMENT_TBB_CONTRAST_ENHANCEMENT_TBB_H_
