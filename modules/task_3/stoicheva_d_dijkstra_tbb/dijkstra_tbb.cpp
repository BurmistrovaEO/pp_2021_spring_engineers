// Copyright 2021 Stoicheva Darya
#include <tbb/tbb.h>
#include <limits>
#include <climits>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "../../../modules/task_3/stoicheva_d_dijkstra_tbb/dijkstra_tbb.h"

void print_graph(const std::vector<int>& graph, const size_t points, const std::string& prefix) {
#ifdef DEBUG_PRINT
    std::cout << "\n" << prefix << "Graph: {\n";
    for (size_t row = 0; row < points; row++) {
        std::cout << prefix << "  ";
        for (size_t column = 0; column < points; column++) {
            std::cout << graph[row * points + column] << ", ";
        }
        std::cout << "\n";
    }
    std::cout << prefix << "}\n";
#endif
}

void print_vector(const std::vector<int>& vector, const size_t size, const std::string& prefix,
    const std::string& label, std::ostream& out) {
    print_vector<int>(vector, size, prefix, label, out);
}


void print_vector(const std::vector<float>& vector, const size_t size, const std::string& prefix,
    const std::string& label, std::ostream& out) {
    print_vector<float>(vector, size, prefix, label, out);
}


void print_vector(const std::vector<double>& vector, const size_t size, const std::string& prefix,
    const std::string& label, std::ostream& out) {
    print_vector<double>(vector, size, prefix, label, out);
}


template<typename T>
void print_vector(const std::vector<T>& vector, const size_t size, const std::string& prefix,
    const std::string& label, std::ostream& out) {
#ifdef DEBUG_PRINT
    std::stringstream buffer;
    buffer << prefix + (label.empty() ? "Vector" : label) + ": {" << std::endl;
    buffer << prefix << "  ";
    for (size_t i = 0; i < size; i++) {
        buffer << vector[i] << ", ";
    }
    buffer << std::endl << prefix << "}" << std::endl;
    out << buffer.str();
#endif
}

int find_unprocessed_point_with_min_distance(const std::vector<int>& graph,
    const std::vector<int>& distances, const std::vector<bool>& processed) {
    int found_point = -1;
    int found_min_distance = MAX_DISTANCE;
    for (size_t point = 0; point < processed.size(); point++) {
        if (!processed[point] && distances[point] < found_min_distance) {
            found_min_distance = distances[point];
            found_point = static_cast<int>(point);
        }
    }
    // std::cout << "FUP: Next point with min distance " << found_min_distance
    //           << " from " << 0 << " is " << found_point << std::endl;
    return found_point;
}

int find_unprocessed_point_with_min_distance_tbb(const std::vector<int>& graph,
    const std::vector<int>& distances, const std::vector<bool>& processed) {

    PointInfo minPoint;
    minPoint.distance = MAX_DISTANCE;
    minPoint.index = -1;
    size_t pointsCount = processed.size();
    size_t pointsPerThread = pointsCount / THREADS_COUNT;

    minPoint = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, static_cast<int>(pointsCount), static_cast<int>(pointsPerThread)),
        minPoint,
        [&](const tbb::blocked_range<int>& range, PointInfo local_minPoint) {
            for (int point = range.begin(); point < range.end(); point++) {
                if (!processed[point] && distances[point] < local_minPoint.distance) {
                    local_minPoint.distance = distances[point];
                    local_minPoint.index = point;
                }
            }
            return local_minPoint;
        },
        [&](PointInfo p1, PointInfo p2) {
            if (p1.distance > p2.distance) {
                return p2;
            }
            return p1;
        });

    return minPoint.index;
}

int process_unprocessed_point(const std::vector<int>& graph,
    std::vector<int>* distances,
    std::vector<bool>* processed, int current_point) {

    int min_distance = MAX_DISTANCE;
    int min_distance_point = -1;
    for (int point = 0; point < static_cast<int>(processed->size()); point++) {
        int start_row_index = current_point * static_cast<int>(processed->size());
        int distance = graph[start_row_index + point];
        if (!(*processed)[point] && distance > 0) {
            int *dp = distances->data() + point;
            int *dcp = distances->data() + current_point;
            *dp =  *dp < *dcp + distance ? *dp : *dcp + distance;  //         std::min(*dp, *dcp + distance);
            if (min_distance > *dp) {
                min_distance = *dp;
                min_distance_point = point;
            }
        }
    }
    (*processed)[current_point] = true;
    // std::cout << "PUP: Next point with min distance " << min_distance << " is "
    //           << min_distance_point + 1 << std::endl;
    return min_distance_point;
}

int process_unprocessed_point_tbb(const std::vector<int>& graph,
    std::vector<int>* distances,
    std::vector<bool>* processed, int current_point) {

    PointInfo minPoint;
    minPoint.distance = MAX_DISTANCE;
    minPoint.index = -1;
    size_t pointsCount = processed->size();
    size_t pointsPerThread = pointsCount / THREADS_COUNT;

    minPoint = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, static_cast<int>(pointsCount), static_cast<int>(pointsPerThread)),
        minPoint,
        [&](const tbb::blocked_range<int>& range, PointInfo local_minPoint) {
            for (int point = range.begin(); point < range.end(); point++) {
                int start_row_index = current_point * static_cast<int>(processed->size());
                int distance = graph[start_row_index + point];
                if (!(*processed)[point] && distance > 0) {
                    int* dp = distances->data() + point;
                    int* dcp = distances->data() + current_point;
                    *dp = *dp < *dcp + distance ? *dp : *dcp + distance;  // std::min(*dp, *dcp + distance);
                    if (local_minPoint.distance > *dp) {
                        local_minPoint.distance = *dp;
                        local_minPoint.index = point;
                    }
                }
            }
            return local_minPoint;
        },
        [&](PointInfo p1, PointInfo p2) {
            if (p1.distance > p2.distance) {
                return p2;
            }
            return p1;
        });

    (*processed)[current_point] = true;
    return minPoint.index;
}

std::vector<int> dijkstra(const std::vector<int>& graph, size_t start, size_t end) {
    if (graph.size() == 0) {
        throw "Error: empty graph";
    }

    size_t points_count = static_cast<size_t>(sqrt(graph.size()));
    if (points_count * points_count != graph.size()) {
        std::cout << "Illegal graph size: expected " << graph.size() << " calculated "
                  << points_count * points_count << std::endl;
        throw "Error: incorrect graph";
    }

    if (start < 1 || end < 1 || start > points_count || end > points_count) {
        std::cout << "Error: illegal start or end" << std::endl;
        throw "Error: illegal start or end";
    }
    start--; end--;

    if (points_count == 1) {
        return { static_cast<int>(start) + 1 };
    }

    constexpr int max_int = MAX_DISTANCE;
    std::vector<int> distances = std::vector<int>(points_count, max_int);
    std::vector<bool> processed = std::vector<bool>(points_count, false);

    distances[start] = 0;

    int next_unprocessed_point = static_cast<int>(start);
    while (next_unprocessed_point >= 0) {
        next_unprocessed_point = process_unprocessed_point(graph,
            &distances, &processed, next_unprocessed_point);
        if (next_unprocessed_point < 0) {
            // std::cout << "- next unprocessed point not given, finding." << std::endl;
            next_unprocessed_point =
                find_unprocessed_point_with_min_distance(graph, distances,
                    processed);
        } else {
            // std::cout << "+ next unprocessed point given, ok." << std::endl;
        }
    }
    print_vector(distances, distances.size(), "distances");

    std::vector<int> path;
    path.push_back(static_cast<int>(end + 1));
    int weight = distances[end];
    int current = static_cast<int>(end);

    while (current != static_cast<int>(start)) {
        for (int i = 0; i < static_cast<int>(points_count); i++) {
            if (graph[current * points_count + i] > 0) {
                int tmp = weight - graph[current * points_count + i];
                if (distances[i] == tmp) {
                    weight = tmp;
                    current = i;
                    path.insert(path.begin(), i + 1);
                    break;
                }
            }
        }
    }

    return path;
}

std::vector<int> dijkstra_tbb(const std::vector<int>& graph, size_t start, size_t end) {
    if (graph.size() == 0)
        throw "Error: empty graph";

    size_t points_count = static_cast<size_t>(sqrt(graph.size()));

    if (points_count * points_count != graph.size()) {
        std::cout << "Illegal graph size: expected " << graph.size() << " calculated "
                  << points_count * points_count << std::endl;
        throw "Error: incorrect graph";
    }

    if (start < 1 || end < 1 || start > points_count || end > points_count) {
        std::cout << "Error: illegal start or end" << std::endl;
        throw "Error: illegal start or end";
    }
    start--; end--;

    if (points_count == 1) {
        return { static_cast<int>(start) + 1 };
    }

    constexpr int max_int = MAX_DISTANCE;
    std::vector<int> distances = std::vector<int>(points_count, max_int);
    std::vector<bool> processed = std::vector<bool>(points_count, false);

    distances[start] = 0;

    tbb::task_scheduler_init init(THREADS_COUNT);
    tbb::mutex mutex;

    int next_unprocessed_point = static_cast<int>(start);
    while (next_unprocessed_point >= 0) {
            next_unprocessed_point = process_unprocessed_point_tbb(graph,
                &distances, &processed, next_unprocessed_point);
            if (next_unprocessed_point < 0) {
                next_unprocessed_point =
                    find_unprocessed_point_with_min_distance_tbb(graph, distances,
                        processed);
            }
    }
    print_vector(distances, distances.size(), "distances:");

    std::vector<int> path;
    path.push_back(static_cast<int>(end + 1));
    int weight = distances[end];
    int current = static_cast<int>(end);

    while (current != static_cast<int>(start)) {
        for (int i = 0; i < static_cast<int>(points_count); i++) {
            if (graph[current * points_count + i] > 0) {
                int tmp = weight - graph[current * points_count + i];
                if (distances[i] == tmp) {
                    weight = tmp;
                    current = i;
                    path.insert(path.begin(), i + 1);
                    break;
                }
            }
        }
    }

    return path;
}
