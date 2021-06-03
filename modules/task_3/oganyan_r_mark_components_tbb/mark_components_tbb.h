// Copyright 2021 Oganyan Robert


#ifndef MODULES_TASK_3_OGANYAN_R_MARK_COMPONENTS_TBB_MARK_COMPONENTS_TBB_H_
#define MODULES_TASK_3_OGANYAN_R_MARK_COMPONENTS_TBB_MARK_COMPONENTS_TBB_H_


#include <utility>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <queue>
#include "../../modules/task_3/oganyan_r_mark_components_tbb/Disjoint_Set_Union.h"




std::pair<std::vector<int>, int> MarkComponentsSeq(std::vector<int> *img,
                                                   int height, int width);

std::pair<std::vector<int>, int> MarkComponentsParTbb(std::vector<int> *img,
                                                   int height, int width, int num_proc);

void bfs(std::vector<int> *img, std::pair<int, int> start,
         int *number, int width, int height);


#endif  //  MODULES_TASK_3_OGANYAN_R_MARK_COMPONENTS_TBB_MARK_COMPONENTS_TBB_H_
