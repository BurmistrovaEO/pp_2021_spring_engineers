// Copyright 2021 Voronin Aleksey
#include "../../../modules/task_4/voronin_a_dijkstra_alg/dijkstra_algorithm.hpp"
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include "../../../3rdparty/unapproved/unapproved.h"

std::vector<int> generateSimpleGraph(int size) {
     std::vector<int> sample(size * size, 2);
     for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                sample[i * size + j] = 1;
            }
        }
     }
     return sample;
}

std::vector<int> dijkstra(const std::vector<int>& graph, int start, int end) {
    if (graph.size() == 0) {
        throw "Empty graph";
    }

    if (start == end)
        return std::vector<int>(1, 0);

    if (start > end) {
        std::swap(start, end);
    }

    int points_count = sqrt(graph.size());

    if (sqrt(graph.size()) != points_count) {
        throw "Wrong size";
    }

    int max_weight = std::numeric_limits<int>::max();
    int min, min_point, tmp;
    std::vector<int> points_len(points_count, max_weight);
    std::vector<int> path;
    std::vector<bool> processed(points_count, false);

    // Align start and end with array indexes
    --start;
    --end;

    points_len[start] = 0;

    do {
        min_point = max_weight;
        min = max_weight;

        // Choose a point to work with
        for (int i = 0; i < points_count; i++) {
            if (!processed[i] && points_len[i] < min) {
                min_point = i;
                min = points_len[i];
            }
        }

        if (min_point != max_weight) {
            for (int i = 0; i < points_count; i++) {
                if (graph[min_point * points_count + i] > 0) {
                    tmp = min + graph[min_point * points_count + i];
                    if (points_len[i] > tmp) {
                        points_len[i] = tmp;
                    }
                }
            }
            processed[min_point] = true;
        }
    } while (min_point < max_weight);

    // Configuring a path
    path.push_back(end + 1);
    int weight = points_len[end];

    while (end != start) {
        for (int i = 0; i < points_count; i++) {
            if (graph[end * points_count + i] < 0) {
                throw "Graph weight can not be less then zero.";
            }
            if (graph[end * points_count + i] > 0) {
                tmp = weight - graph[end * points_count + i];
                if (points_len[i] == tmp) {
                    weight = tmp;
                    end = i;
                    path.push_back(i + 1);
                }
            }
        }
    }

    return path;
}


std::vector<int> dijkstraParallel(const std::vector<int>& graph, int start, int end) {
    if (graph.size() == 0) {
        throw "Empty graph";
    }

    if (start == end)
        return std::vector<int>(1, 0);

    if (start > end) {
        std::swap(start, end);
    }

    int points_count = sqrt(graph.size());

    if (sqrt(graph.size()) != points_count) {
        throw "Wrong size";
    }

    int max_weight = std::numeric_limits<int>::max();
    int min, min_point, tmp;
    std::vector<int> points_len(points_count, max_weight);
    std::vector<int> path;
    std::vector<bool> processed(points_count, false);

    // Align start and end with array indexes
    --start;
    --end;

    points_len[start] = 0;

    do {
        min_point = max_weight;
        min = max_weight;

        // Choose a point to work with
        for (int i = 0; i < points_count; i++) {
            if (!processed[i] && points_len[i] < min) {
                min_point = i;
                min = points_len[i];
            }
        }

        std::vector<std::thread> threads;
        int thread_count = 4;
        if (min_point != max_weight) {
            for (int k = 0; k < thread_count ; k++) {
            threads.emplace_back([&graph, &points_count, &min, &points_len, &min_point, thread_count, k]() {
                int temp;
                for (int i = points_count / thread_count * k; i < points_count/thread_count * (k + 1); i++) {
                    if (graph[min_point * points_count + i] > 0) {
                        temp = min + graph[min_point * points_count + i];
                        if (points_len[i] > temp) {
                            points_len[i] = temp;
                        }
                    }
                }
            });
            }
            for (auto& t : threads) {
            t.join();
            }
            processed[min_point] = true;
        }
    } while (min_point < max_weight);

    // Configuring a path
    path.push_back(end + 1);
    int weight = points_len[end];

    while (end != start) {
        for (int i = 0; i < points_count; i++) {
            if (graph[end * points_count + i] < 0) {
                throw "Graph weight can not be less then zero.";
            }
            if (graph[end * points_count + i] > 0) {
                tmp = weight - graph[end * points_count + i];
                if (points_len[i] == tmp) {
                    weight = tmp;
                    end = i;
                    path.push_back(i + 1);
                }
            }
        }
    }

    return path;
}
