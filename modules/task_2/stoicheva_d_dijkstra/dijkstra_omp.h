// Copyright 2021 Stoicheva Darya

#ifndef MODULES_TASK_2_STOICHEVA_D_DIJKSTRA_DIJKSTRA_OMP_H_
#define MODULES_TASK_2_STOICHEVA_D_DIJKSTRA_DIJKSTRA_OMP_H_

#include <vector>
#include <string>
#include <random>
#include <limits>

// #define DEBUG_PRINT

void print_graph(const std::vector<int>& graph, const size_t points, const std::string& prefix = {});
void print_vector(const std::vector<int>& vector, const size_t size, const std::string& prefix = {},
    const std::string& label = {}, std::ostream& out = std::cout);

template<typename T>
void print_vector(const std::vector<T>& vector, const size_t size, const std::string& prefix = {},
    const std::string& label = {}, std::ostream& out = std::cout);

std::vector<int> dijkstra(const std::vector<int>& graph, size_t start, size_t end);
std::vector<int> dijkstra_omp(const std::vector<int>& graph, size_t start, size_t end);

#endif  // MODULES_TASK_2_STOICHEVA_D_DIJKSTRA_DIJKSTRA_OMP_H_
