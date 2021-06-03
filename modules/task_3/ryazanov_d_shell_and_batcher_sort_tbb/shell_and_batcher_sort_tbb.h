// Copyright 2021 Dmitry Ryazanov
#ifndef MODULES_TASK_3_RYAZANOV_D_SHELL_AND_BATCHER_SORT_TBB_SHELL_AND_BATCHER_SORT_TBB_H_
#define MODULES_TASK_3_RYAZANOV_D_SHELL_AND_BATCHER_SORT_TBB_SHELL_AND_BATCHER_SORT_TBB_H_

#include <tbb/tbb.h>
#include <algorithm>
#include <stdexcept>
#include <vector>

std::vector<int> getRandomVector(int sz);
std::vector<int> SortShell(const std::vector<int>& vec);
std::vector<int> SortShellTBB(const std::vector<int>& vec);
void print(std::vector<int> const& a);
std::vector<int> sh(const std::vector<int>& vec, int l, int r);
std::vector<int> uns(const std::vector<int>& vec, int l, int r);
std::vector<int> unsTBB(const std::vector<int>& vec, int l, int r);
std::vector<int> BatcherSort(const std::vector<int>& vec);
std::vector<int> BatcherSortTBB(const std::vector<int>& vec);

#endif  // MODULES_TASK_3_RYAZANOV_D_SHELL_AND_BATCHER_SORT_TBB_SHELL_AND_BATCHER_SORT_TBB_H_
