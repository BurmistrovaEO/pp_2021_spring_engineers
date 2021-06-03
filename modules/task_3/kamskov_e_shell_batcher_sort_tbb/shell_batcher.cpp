// Copyright 2021 Kamskov Eugene

#include "../../../modules/task_3/kamskov_e_shell_batcher_sort_tbb/shell_batcher.h"
#include <tbb/tbb.h>
#include <omp.h>

std::vector<int> getRandomVector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> res(size);
    for (int i = 0; i < size; ++i)
        res[i] = gen() % 100;
    return res;
}

void print(std::vector<int> const& a) {
    for (size_t i = 0; i < a.size(); ++i) std::cout << a[i] << ' ';
    std::cout << std::endl;
}

std::vector<std::vector<int>> splitV(const std::vector<int>& vec, size_t n) {
    size_t del = vec.size() / n;
    size_t bal = vec.size() % n;
    size_t begin = 0, end = 0;
    std::vector<std::vector<int>> res;
    for (size_t i = 0; i < n; ++i) {
        if (bal > 0) {
            end = end + (del + !!(bal--));
        } else {
            end = end + del;
        }
        res.push_back(std::vector<int>(vec.begin() + begin, vec.begin() + end));
        begin = end;
    }
  return res;
}

std::vector<int> shellSort(const std::vector<int>& v, int size) {
    int del, i, j;
    std::vector<int>res(v);
    for (del = size / 2; del > 0; del = del / 2) {
        for (i = del; i < static_cast<int>(res.size()); i++) {
            for (j = i - del; j >= 0 && res[j] > res[j + del]; j = j - del) {
                std::swap(res[j], res[j + del]);
            }
        }
    }
    return res;
}

std::vector<int> batchShellSort(const std::vector<int>& v, int size) {
    std::vector<int>fHalf;
    std::copy(v.begin(), v.end() - size / 2, inserter(fHalf, fHalf.begin()));
    fHalf = shellSort(fHalf, fHalf.size());
    std::vector<int> sHalf;
    std::copy(v.end() - size / 2, v.end(), inserter(sHalf, sHalf.begin()));
    sHalf = shellSort(sHalf, sHalf.size());

    std::vector<int>even = evenBatch(fHalf, sHalf), odd = oddBatch(fHalf, sHalf);
    std::vector<int>res = mergeBatch(even, odd);
    return res;
}

std::vector<int> batchShellSortTbb(const std::vector<int>& A, const int n, int size) {
    std::vector<std::vector<int>> vec = splitV(A, n);
    tbb::task_scheduler_init init(n);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, vec.size(), 1),
        [&vec](const tbb::blocked_range<size_t>& r) {
        int begin = r.begin(), end = r.end();
        for (int i = begin; i != end; ++i)
            vec[i] = shellSort(vec[i], vec[i].size());
        }, tbb::simple_partitioner());

    init.terminate();
    std::vector<int> res = vec[0];
    std::vector<int> odd, even;

    for (int i = 1; i < n; i++) {
        odd = oddBatch(res, vec[i]);
        even = evenBatch(res, vec[i]);
        res = mergeBatch(even, odd);
    }
    return res;
}

std::vector<int> evenBatch(const std::vector<int>& A, const std::vector<int>& B) {
    std::vector <int> res(A.size() / 2 + B.size() / 2 + A.size() % 2 + B.size() % 2);
    size_t a = 0, b = 0, i = 0;
    while ((a < A.size()) && (b < B.size())) {
        if (A[a] <= B[b]) {
            res[i] = A[a];
            a = a + 2;
        } else {
            res[i] = B[b];
            b = b + 2;
        }
        ++i;
    }
    if (a >= A.size()) {
        for (size_t j = b; j < B.size(); j = j + 2) {
            res[i] = B[j];
            ++i;
        }
    } else {
        for (size_t j = a; j < A.size(); j = j + 2) {
            res[i] = A[j];
            ++i;
        }
    }
    return res;
}

std::vector<int> oddBatch(const std::vector<int>& A, const std::vector<int>& B) {
    std::vector <int> res(A.size() / 2 + B.size() / 2);
    size_t a = 1, b = 1, i = 0;
    while ((a < A.size()) && (b < B.size())) {
        if (A[a] <= B[b]) {
            res[i] = A[a];
            a = a + 2;
        } else {
            res[i] = B[b];
            b = b + 2;
        }
        ++i;
    }
    if (a >= A.size()) {
        for (size_t j = b; j < B.size(); j = j + 2) {
            res[i] = B[j];
            ++i;
        }
    } else {
        for (size_t j = a; j < A.size(); j = j + 2) {
            res[i] = A[j];
            ++i;
        }
    }
    return res;
}

std::vector<int> mergeBatch(const std::vector<int>& even, const std::vector<int>& odd) {
    size_t sizeE = even.size(), sizeO = odd.size(), j, size = sizeE + sizeO, i = 0, ev = 0, od = 0;
    int tmp;
    std::vector<int> res(size);
    while ((ev < sizeE) && (od < sizeO)) {
        res[i] = even[ev];
        res[i + 1] = odd[od];
        i = i + 2;
        ++ev;
        ++od;
    }
    if ((od >= sizeO) && (ev < sizeE)) {
        for (j = i; j < size; j++) {
            res[j] = even[ev];
            ++ev;
        }
    }
    for (size_t i = 0; i < size - 1; i++) {
        if (res[i] > res[i + 1]) {
            tmp = res[i];
            res[i] = res[i + 1];
            res[i + 1] = tmp;
        }
    }
    return res;
}

