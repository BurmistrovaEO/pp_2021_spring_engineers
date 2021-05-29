// Copyright 2021 Stoicheva Darya
#include <omp.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "../../../modules/task_2/stoicheva_d_dijkstra/dijkstra_omp.h"

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
    int found_min_distance = std::numeric_limits<int>::max();
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

int find_unprocessed_point_with_min_distance_omp(const std::vector<int>& graph,
    const std::vector<int>& distances, const std::vector<bool>& processed) {
    int found_point = -1;
    int found_min_distance = std::numeric_limits<int>::max();
    #pragma omp parallel
    {
        int local_found_point = found_point;
        int local_found_min_distance = found_min_distance;
        #pragma omp for
        for (int point = 0; point < static_cast<int>(processed.size()); point++) {
            if (!processed[point] && distances[point] < local_found_min_distance) {
                local_found_min_distance = distances[point];
                local_found_point = point;
            }
        }
        #pragma omp critical
        {
            if (found_min_distance > local_found_min_distance) {
                found_min_distance = local_found_min_distance;
                found_point = local_found_point;
            }
        }
    }
    return found_point;
}

int process_unprocessed_point(const std::vector<int>& graph,
    std::vector<int>* distances,
    std::vector<bool>* processed, int current_point) {

    int min_distance = std::numeric_limits<int>::max();
    int min_distance_point = -1;
    for (int point = 0; point < static_cast<int>(processed->size()); point++) {
        int start_row_index = current_point * static_cast<int>(processed->size());
        int distance = graph[start_row_index + point];
        if (!(*processed)[point] && distance > 0) {
            int *dp = distances->data() + point;
            int *dcp = distances->data() + current_point;
            *dp = std::min(*dp, *dcp + distance);
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

int process_unprocessed_point_omp(const std::vector<int>& graph,
    std::vector<int>* distances,
    std::vector<bool>* processed, int current_point) {

    int min_distance = std::numeric_limits<int>::max();
    int min_distance_point = -1;
    #pragma omp parallel
    {
        int local_min_distance = std::numeric_limits<int>::max();
        int local_min_distance_point = -1;
        #pragma omp for
        for (int point = 0; point < static_cast<int>(processed->size()); point++) {
            int start_row_index = current_point * static_cast<int>(processed->size());
            int distance = graph[start_row_index + point];
            if (!(*processed)[point] && distance > 0) {
                int* dp = distances->data() + point;
                int* dcp = distances->data() + current_point;
                *dp = std::min(*dp, *dcp + distance);
                if (local_min_distance > *dp) {
                    local_min_distance = *dp;
                    local_min_distance_point = point;
                }
            }
        }
        #pragma omp critical
        {
            if (min_distance > local_min_distance) {
                min_distance = local_min_distance;
                min_distance_point = local_min_distance_point;
            }
        }
    }
    (*processed)[current_point] = true;
    return min_distance_point;
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

    constexpr int max_int = std::numeric_limits<int>::max();
    std::vector<int> distances = std::vector<int>(points_count, max_int);
    std::vector<bool> processed = std::vector<bool>(points_count, false);

    distances[start] = 0;

    int next_unprocessed_point = start;
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
    path.push_back(end + 1);
    int weight = distances[end];
    int current = end;

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

std::vector<int> dijkstra_omp(const std::vector<int>& graph, size_t start, size_t end) {
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

    constexpr int max_int = std::numeric_limits<int>::max();
    std::vector<int> distances = std::vector<int>(points_count, max_int);
    std::vector<bool> processed = std::vector<bool>(points_count, false);

    distances[start] = 0;

    int next_unprocessed_point = start;
    while (next_unprocessed_point >= 0) {
            next_unprocessed_point = process_unprocessed_point_omp(graph,
                &distances, &processed, next_unprocessed_point);
            if (next_unprocessed_point < 0) {
                next_unprocessed_point =
                    find_unprocessed_point_with_min_distance_omp(graph, distances,
                        processed);
            }
    }
    print_vector(distances, distances.size(), "distances:");

    std::vector<int> path;
    path.push_back(end + 1);
    int distance = distances[end];
    int current = end;
    bool found = false;

    while (current != static_cast<int>(start)) {
        found = false;

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(points_count); i++) {
            // if (found) {
            //     std::cout << "################ Found and not skipped #################" << std::endl;
            // }
            if (!found && graph[current * points_count + i] > 0) {
                int tmp_distance = distance - graph[current * points_count + i];
                if (distances[i] == tmp_distance) {
                    #pragma omp critical
                    {
                        distance = tmp_distance;
                        current = i;
                        path.insert(path.begin(), i + 1);
                        found = true;
                        // std::cout << "################ Found - now request for skip #################" << std::endl;
                    }
                }
            }
        }
    }

    return path;
}
