// Copyright 2021 Oganyan Robert

#ifndef MODULES_TASK_4_OGANYAN_R_MARK_COMPONENTS_STD_DISJOINT_SET_UNION_H_
#define MODULES_TASK_4_OGANYAN_R_MARK_COMPONENTS_STD_DISJOINT_SET_UNION_H_

#include <vector>
#include <utility>
#include <algorithm>

template<typename T>
class Disjoint_Set_Union {
 private:
    std::vector<T> rank;
    std::vector<T> parent;

 public:
    explicit Disjoint_Set_Union(int size) : rank(size), parent(size) {
    }

    void make_set(int vertex) {
        parent[vertex] = vertex;
        rank[vertex] = 0;
    }


    void Init() {
        for (std::size_t vertex = 0; vertex < rank.size(); ++vertex) {
            make_set(vertex);
        }
    }


    void Init(int size) {
        rank.resize(size);
        parent.resize(size);
        for (int vertex = 0; vertex < size; ++vertex) {
            make_set(vertex);
        }
    }

    int find_set(int vertex) {
        //  Bug here (parent[vertex] = [parent[paren[vertex]])
        if (vertex == parent[vertex]) {
            return vertex;
        }

        return parent[vertex] = find_set(parent[vertex]);
    }

    int find_set(int vertex, int last) {
        if (vertex == parent[vertex]) {
            return vertex;
        }

        if (last == parent[vertex]) {
            return std::min(last, vertex);
        }

        return parent[vertex] = find_set(parent[vertex], vertex);
    }

    void union_sets(int fi_union, int se_union) {
        fi_union = find_set(fi_union, fi_union);
        se_union = find_set(se_union, se_union);
        if (fi_union != se_union) {
            if (rank[fi_union] < rank[se_union])
                std::swap(fi_union, se_union);
            parent[se_union] = fi_union;
            if (rank[fi_union] == rank[se_union])
                ++rank[fi_union];
        }
    }

    const std::vector<T> get_rank() {
        return this->rank;
    }

    const std::vector<T> get_parent() {
        return this->parent;
    }
};

#endif  //  MODULES_TASK_4_OGANYAN_R_MARK_COMPONENTS_STD_DISJOINT_SET_UNION_H_
