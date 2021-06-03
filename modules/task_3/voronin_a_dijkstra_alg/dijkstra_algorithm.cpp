// Copyright 2021 Voronin Aleksey
#include "../../../modules/task_3/voronin_a_dijkstra_alg/dijkstra_algorithm.hpp"
#include <tbb/tbb.h>
#include <vector>
#include <limits>
#include <cmath>
#define THREADS 4

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

    int max_weight =  INT32_MAX;
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

    int max_weight =  INT32_MAX;
    int min, min_point, tmp;
    std::vector<int> points_len(points_count, max_weight);
    std::vector<int> path;
    std::vector<bool> processed(points_count, false);
    tbb::mutex mutex;
    tbb::task_scheduler_init(THREADS);
    std::vector<int> min_vals(2, INT32_MAX);

    --start;
    --end;

    points_len[start] = 0;

    do {
        min_vals = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, points_count),
            std::vector<int>(2) = {INT32_MAX, INT32_MAX},
            [&](const tbb::blocked_range<int>& v, std::vector<int> local_min_vals) {
                for (int i = v.begin(); i < v.end(); i++) {
                    if (!processed[i] && points_len[i] < local_min_vals[0]) {
                        local_min_vals[0] = points_len[i];
                        local_min_vals[1] = i;
                    }
                }
                return local_min_vals;
            },
            [&](std::vector<int> x, std::vector<int> y) {
                if (x[0] < y[0]) {
                    return x;
                }
                return y;
            });

        min_point = min_vals[1];
        min = min_vals[0];

        if (min_point != max_weight) {
            tbb::parallel_for(
                tbb::blocked_range<int>(0, points_count),
                [&](const tbb::blocked_range<int>& v) {
                    for (int i = v.begin(); i < v.end(); i++) {
                        if (graph[min_point * points_count + i] > 0) {
                            mutex.lock();
                            tmp = min + graph[min_point * points_count + i];
                            if (points_len[i] > tmp) {
                                points_len[i] = tmp;
                            }
                            mutex.unlock();
                        }
                    }
                });
            processed[min_point] = true;
        }
    } while (min_point < max_weight);

    path.push_back(end + 1);
    int weight = points_len[end];

    while (end != start) {
        tbb::parallel_for(
            tbb::blocked_range<int>(0, points_count),
            [&](const tbb::blocked_range<int>& v) {
                for (int i = v.begin(); i < v.end(); i++) {
                    if (graph[end * points_count + i] > 0) {
                        tmp = weight - graph[end * points_count + i];
                        if (points_len[i] == tmp) {
                            weight = tmp;
                            end = i;
                            mutex.lock();
                            path.push_back(i + 1);
                            mutex.unlock();
                        }
                    }
                }
            });
    }

    return path;
}
