// Copyright 2020 Stoicheva Darya
#include <vector>
#include <string>
#include <utility>
#include <random>
#include <iostream>
#include <sstream>
#include "../../../modules/task_4/stoicheva_d_dijkstra_std/stoicheva_d_dijkstra_std.h"
#include "../../../3rdparty/unapproved/unapproved.h"

std::vector<int> getRandomVector(int sz) {
    std::random_device dev;
    std::mt19937 gen(dev());
    std::vector<int> vec(sz);
    for (int  i = 0; i < sz; i++) { vec[i] = gen() % 100; }
    return vec;
}

std::mutex my_mutex;

void atomOps(std::vector<int> vec, const std::string& ops, std::promise<int> &&pr) {
    const int  sz = vec.size();
    int reduction_elem = 0;
    if (ops == "+") {
        for (int  i = 0; i < sz; i++) {
            std::lock_guard<std::mutex> my_lock(my_mutex);
            reduction_elem += vec[i];
        }
    } else if (ops == "-") {
        for (int  i = 0; i < sz; i++) {
            std::lock_guard<std::mutex> my_lock(my_mutex);
            reduction_elem -= vec[i];
        }
    }
    pr.set_value(reduction_elem);
}



int getParallelOperations(std::vector<int> vec, const std::string& ops) {
    const int nthreads = std::thread::hardware_concurrency();
    const int delta = (vec.end() - vec.begin()) / nthreads;

    std::promise<int> *promises = new std::promise<int>[nthreads];
    std::future<int> *futures = new std::future<int>[nthreads];
    std::thread *threads = new std::thread[nthreads];

    int reduction_elem = 0;
    for (int i = 0; i < nthreads; i++) {
        futures[i] = promises[i].get_future();
        std::vector<int> tmp_vec(
            vec.begin() + i * delta,
            vec.begin() + (i + 1) * delta);
        threads[i] = std::thread(atomOps, tmp_vec, ops, std::move(promises[i]));
        threads[i].join();
        reduction_elem += futures[i].get();
    }

    delete []promises;
    delete []futures;
    delete []threads;
    return reduction_elem;
}

int getSequentialOperations(std::vector<int> vec, const std::string& ops) {
    const int  sz = vec.size();
    int reduction_elem = 0;
    if (ops == "+") {
        for (int  i = 0; i < sz; i++) {
            reduction_elem += vec[i];
        }
    } else if (ops == "-") {
        for (int  i = 0; i < sz; i++) {
            reduction_elem -= vec[i];
        }
    }
    return reduction_elem;
}

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


void atom_find_unprocessed_point_with_min_distance_std(const std::vector<int>& graph, int start, int end,
    const std::vector<int>& distances, const std::vector<bool>& processed, std::promise<PointInfo> &&pr) {
    const int size = graph.size();
    PointInfo local_minPoint;
    local_minPoint.distance = MAX_DISTANCE;
    local_minPoint.index = -1;
    for (int point = start; point < end; point++) {
        if (!processed[point] && distances[point] < local_minPoint.distance) {
            local_minPoint.distance = distances[point];
            local_minPoint.index = point;
        }
    }
    pr.set_value(local_minPoint);
}


int find_unprocessed_point_with_min_distance_std(const std::vector<int>& graph,
    const std::vector<int>& distances, const std::vector<bool>& processed) {


    size_t pointsCount = processed.size();


    const int nthreads = std::thread::hardware_concurrency();  //  THREADS_COUNT
    size_t pointsPerThread = pointsCount / nthreads;
    const int delta = (graph.end() - graph.begin()) / nthreads;
    int last = 0;
    if (pointsCount > nthreads * delta) {
        last = pointsCount - nthreads * delta;
    }

    std::promise<PointInfo> *promises = new std::promise<PointInfo>[nthreads];
    std::future<PointInfo> *futures = new std::future<PointInfo>[nthreads];
    std::thread *threads = new std::thread[nthreads];

#ifdef DEBUG_PRINT
    std::cout << "Threads: " << nthreads << ", Points per thread: " << pointsPerThread
        << ", last: " << last << std::endl;
#endif


    PointInfo minPoint;
    minPoint.distance = MAX_DISTANCE;
    minPoint.index = -1;
    for (int i = 0; i < nthreads; i++) {
        futures[i] = promises[i].get_future();
        int start = i * pointsPerThread;
        int end = i < nthreads - 1 ? (i + 1) * pointsPerThread : pointsCount;
#ifdef DEBUG_PRINT
        std::cout << "Thread: " << i << ", start: " << start << ", end: " << end << std::endl;
#endif
        threads[i] = std::thread(atom_find_unprocessed_point_with_min_distance_std,
            std::ref(graph), start, end, std::ref(distances), std::ref(processed), std::move(promises[i]));
        threads[i].detach();
    }

    for (int i = 0; i < nthreads; i++) {
        PointInfo threadMinPoint = futures[i].get();
        if (minPoint.distance > threadMinPoint.distance) {
            minPoint.distance = threadMinPoint.distance;
            minPoint.index = threadMinPoint.index;
        }
    }

    delete []promises;
    delete []futures;
    delete []threads;

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
            if (*dp > *dcp + distance) {
                *dp = *dcp + distance;
            }
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

void atom_process_unprocessed_point_std(
    const std::vector<int>& graph, const int start, const int end,
    std::vector<int>* distances, std::vector<bool>* processed,
    const int current_point, std::promise<PointInfo>&& pr) {
    int start_row_index = current_point * static_cast<int>(processed->size());
    int* dcp = &(*distances)[current_point];
    PointInfo local_minPoint;
    local_minPoint.distance = MAX_DISTANCE;
    local_minPoint.index = -1;
    for (int point = start; point < end; point++) {
        int distance = graph[start_row_index + point];
        if (!(*processed)[point] && distance > 0) {
            int* dp = &(*distances)[point];
            std::lock_guard<std::mutex> my_lock(my_mutex);
            *dp = *dp < *dcp + distance ? *dp : *dcp + distance;  // std::min(*dp, *dcp + distance);
            if (local_minPoint.distance > *dp) {
                local_minPoint.distance = *dp;
                local_minPoint.index = point;
            }
        }
    }
    pr.set_value(local_minPoint);
}

int process_unprocessed_point_std(const std::vector<int>& graph,
    std::vector<int>* distances,
    std::vector<bool>* processed, int current_point) {

    const int nthreads = std::thread::hardware_concurrency();  //  THREADS_COUNT
    const int delta = (graph.end() - graph.begin()) / nthreads;

    size_t pointsCount = processed->size();
    size_t pointsPerThread = pointsCount / nthreads;

    std::promise<PointInfo> *promises = new std::promise<PointInfo>[nthreads];
    std::future<PointInfo> *futures = new std::future<PointInfo>[nthreads];
    std::thread *threads = new std::thread[nthreads];

    PointInfo minPoint;
    minPoint.distance = MAX_DISTANCE;
    minPoint.index = -1;
    for (int i = 0; i < nthreads; i++) {
        futures[i] = promises[i].get_future();
        int start = i * pointsPerThread;
        int end = i < nthreads - 1 ? (i + 1) * pointsPerThread : pointsCount;
        threads[i] = std::thread(atom_process_unprocessed_point_std, std::cref(graph), start, end, std::ref(distances),
            std::ref(processed), current_point, std::move(promises[i]));
        threads[i].detach();
    }

    for (int i = 0; i < nthreads; i++) {
        PointInfo threadMinPoint = futures[i].get();
        if (minPoint.distance > threadMinPoint.distance) {
            minPoint.distance = threadMinPoint.distance;
            minPoint.index = threadMinPoint.index;
        }
    }

    delete []promises;
    delete []futures;
    delete []threads;

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
            //  std::cout << "- next unprocessed point not given, finding." << std::endl;
            next_unprocessed_point =
                find_unprocessed_point_with_min_distance(graph, distances,
                    processed);
        } else {
            //  std::cout << "+ next unprocessed point given, ok." << std::endl;
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

std::vector<int> dijkstra_std(const std::vector<int>& graph, size_t start, size_t end) {
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
    int next_unprocessed_point = static_cast<int>(start);
    while (next_unprocessed_point >= 0) {
            next_unprocessed_point = process_unprocessed_point_std(graph,
                &distances, &processed, next_unprocessed_point);
            if (next_unprocessed_point < 0) {
                //  std::cout << "- next unprocessed point not given, finding." << std::endl;
                next_unprocessed_point =
                    find_unprocessed_point_with_min_distance_std(graph, distances,
                        processed);
            } else {
                //  std::cout << "+ next unprocessed point given, ok." << std::endl;
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
