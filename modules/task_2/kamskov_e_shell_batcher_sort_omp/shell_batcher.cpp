// Copyright 2021 Kamskov Eugene
#include "../../../modules/task_2/kamskov_e_shell_batcher_sort_omp/shell_batcher.h"

std::vector<int> getRandomVector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(size);
    for (int i = 0; i < size; ++i)
        vec[i] = gen() % 100;
    return vec;
}
std::vector<std::vector<int>> split(const std::vector<int>& vec, size_t n) {
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

std::vector<int> batchShellSortOmp(const std::vector<int>& v, const int n, int size) {
    std::vector<std::vector<int>> vec = split(v, n);
    std::vector<int> res;
    #pragma omp parallel shared(vec)
    {
    #pragma omp for
        for (int i = 0; i < static_cast<int> (vec.size()); i++) {
            vec[i] = shellSort(vec[i], vec[i].size());
        }
    }
    res = mergeOmp(vec, n, size);
    return res;
}

std::vector<int> mergeOmp(const std::vector<std::vector<int>>& A, const int n, int size) {
    std::vector<int> res = A[0], odd, even;

    for (int i = 1; i < n; i++) {
#pragma omp parallel sections num_threads(1)
        {
#pragma omp section
            {
                odd = oddBatch(res, A[i]);
            }
#pragma omp section
            {
                even = evenBatch(res, A[i]);
            }
#pragma omp section
            {
                res = mergeBatch(even, odd);
            }
        }
    }
    return res;
}

std::vector<int> evenBatch(const std::vector<int>& A, const std::vector<int>& B) {
    int size_res = A.size() / 2 + B.size() / 2 + A.size() % 2 + B.size() % 2;
    std::vector <int> res(size_res);
    size_t iA = 0, iB = 0, i = 0;
    while ((iA < A.size()) && (iB < B.size())) {
        if (A[iA] <= B[iB]) {
            res[i] = A[iA];
            iA += 2;
        } else {
            res[i] = B[iB];
            iB += 2;
        }
        ++i;
    }
    if (iA >= A.size()) {
        for (size_t j = iB; j < B.size(); j = j + 2) {
            res[i] = B[j];
            ++i;
        }
    } else {
        for (size_t j = iA; j < A.size(); j = j + 2) {
            res[i] = A[j];
            ++i;
        }
    }
    return res;
}

std::vector<int> oddBatch(const std::vector<int>& A, const std::vector<int>& B) {
    int size_res = A.size() / 2 + B.size() / 2;
    std::vector <int> res(size_res);
    size_t iA = 1, iB = 1, i = 0;
    while ((iA < A.size()) && (iB < B.size())) {
        if (A[iA] <= B[iB]) {
            res[i] = A[iA];
            iA = iA + 2;
        } else {
            res[i] = B[iB];
            iB = iB + 2;
        }
        ++i;
    }
    if (iA >= A.size()) {
        for (size_t j = iB; j < B.size(); j = j + 2) {
            res[i] = B[j];
            ++i;
        }
    } else {
        for (size_t j = iA; j < A.size(); j = j + 2) {
            res[i] = A[j];
            ++i;
        }
    }
    return res;
}

std::vector<int> mergeBatch(const std::vector<int>& even, const std::vector<int>& odd) {
    size_t sizeE = even.size();
    size_t sizeO = odd.size();
    size_t size = sizeE + sizeO;
    std::vector<int> res(size);
    size_t i = 0, e = 0, o = 0;
    while ((e < sizeE) && (o < sizeO)) {
        res[i] = even[e];
        res[i + 1] = odd[o];
        i = i + 2;
        ++e;
        ++o;
    }
    if ((o >= sizeO) && (e < sizeE)) {
        for (size_t l = i; l < size; l++) {
            res[l] = even[e];
            ++e;
        }
    }
    for (size_t i = 0; i < size - 1; i++) {
        if (res[i] > res[i + 1]) {
            std::swap(res[i], res[i + 1]);
        }
    }
    return res;
}

void print(std::vector<int> const& a) {
    for (size_t i = 0; i < a.size(); ++i) std::cout << a[i] << ' ';
    std::cout << std::endl;
}
