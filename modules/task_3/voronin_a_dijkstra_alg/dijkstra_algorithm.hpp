// Copyright 2021 Voronin Aleksey

#ifndef MODULES_TASK_3_VORONIN_A_DIJKSTRA_ALG_DIJKSTRA_ALGORITHM_HPP_
#define MODULES_TASK_3_VORONIN_A_DIJKSTRA_ALG_DIJKSTRA_ALGORITHM_HPP_
#include <vector>

std::vector<int> generateSimpleGraph(int size);
std::vector<int> dijkstra(const std::vector<int>& graph, int start, int end);
std::vector<int> dijkstraParallel(const std::vector<int>& graph, int start, int end);

#endif  // MODULES_TASK_3_VORONIN_A_DIJKSTRA_ALG_DIJKSTRA_ALGORITHM_HPP_
