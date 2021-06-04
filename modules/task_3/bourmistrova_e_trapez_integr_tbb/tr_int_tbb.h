// Copyright 2021 Ekaterina Burmistrova
#ifndef MODULES_TASK_3_BOURMISTROVA_E_TRAPEZ_INTEGR_TBB_TR_INT_TBB_H_
#define MODULES_TASK_3_BOURMISTROVA_E_TRAPEZ_INTEGR_TBB_TR_INT_TBB_H_
#include <vector>
#include <string>
#include <functional>
#include <utility>

// std::vector<int> getRandomVector(int  sz);
double SolveParallel(const std::vector<std::pair<int, int>>& bord,
	std::function<double(double, double, double)> f);
double SolveParallelSum(const std::vector<std::pair<int, int>>& bord,
	std::function<double(double, double, double)> f);

#endif  // MODULES_TASK_3_BOURMISTROVA_E_TRAPEZ_INTEGR_TBB_TR_INT_TBB_H_


