// Copyright 2018 Nesterov Alexander
#include <tbb/tbb.h>
#include <vector>
#include <string>
#include <random>
#include <functional>
#include <numeric>
#include "../../../modules/task_3/ptitsyn_p_bbb_sort_tbb/bbb_sort.h"

double RandomDouble::Next() {
    static RandomDouble rand = RandomDouble();
    return rand.dist_(rand.gen_);
}

RandomDouble::RandomDouble() : gen_(std::random_device()()), dist_(0, 0e3) {}

std::vector<double> random_double_array(size_t size) {
    std::vector<double> arr(size);

    for (size_t i = 0; i < size; ++i) {
        arr[i] = RandomDouble::Next();
    }

    return arr;
}

std::vector<double> transit(std::vector<double>
    inputNumbers, std::vector<double> loc, const int size, int add) {
    unsigned char* pmem = reinterpret_cast<unsigned char*>(inputNumbers.data());
    int counters[256] = { 0 };
    int sum = 0;

    for (int i = 0; i < size; i++) {
        counters[pmem[8 * i + add]]++;
    }
    for (int i = 0; i < 256; i++) {
        int tmp = counters[i];
        counters[i] = sum;
        sum += tmp;
    }
    for (int i = 0; i < size; i++) {
        int index = 8 * i + add;
        loc[counters[pmem[index]]] = inputNumbers[i];
        counters[pmem[index]]++;
    }
    return loc;
}

std::vector<double> loc_sort(std::vector<double> inputNumbers, const int size) {
    std::vector<double> loc(size, 0);
    std::vector<double> tmp(size, 0);

    for (int i = 0; i < 8; i++) {
        loc = transit(inputNumbers, loc, size, i);
        tmp = inputNumbers;
        inputNumbers = loc;
        loc = tmp;
    }
    return inputNumbers;
}

std::vector<double> BatcherSort(const std::vector<double>& vec) {
    int size = vec.size();

    if (vec.size() == 2) loc_sort(vec, size);
    if (vec.size() == 1) return vec;
    std::vector<double> vecAuxLeft, vecAuxRight, res;
    res = vec;
    vecAuxLeft.reserve(vec.size() / 2);
    vecAuxRight.reserve(vec.size() / 2);
    res = unsTBB(res, 0, size);
    vecAuxLeft.insert(vecAuxLeft.end(), res.begin(),
        res.begin() + res.size() / 2);
    vecAuxRight.insert(vecAuxRight.end(), res.begin() + res.size() / 2,
        res.end());
    vecAuxLeft = BatcherSort(vecAuxLeft);
    vecAuxRight = BatcherSort(vecAuxRight);
    res.erase(res.begin(), res.end());
    res.shrink_to_fit();
    res.insert(res.end(), vecAuxLeft.begin(), vecAuxLeft.end());
    res.insert(res.end(), vecAuxRight.begin(), vecAuxRight.end());
    res = sh(res, 0, size);
    return res;
}

std::vector<double> uns(const std::vector<double>& vec, int l, int r) {
    std::vector<double> auxLeft, auxRight, res;
    res = vec;
    if (res.size() % 2 != 0) {
        int lValue = res.at(res.size() - 1), lPos = res.size() - 1;

        res.erase(res.end() - 1);

        for (int j = l; j < r - 1; j += 2) {
            auxLeft.push_back(res.at(j));
            auxRight.push_back(res.at(j + 1));
        }

        if (lPos % 2 == 0)
            auxLeft.push_back(lValue);
        else
            auxRight.push_back(lValue);
    } else {
        for (int j = l; j < r; j += 2) {
            auxLeft.push_back(res.at(j));
            auxRight.push_back(res.at(j + 1));
        }
    }

    res.erase(res.begin(), res.end());
    res.shrink_to_fit();

    res.reserve(auxLeft.size() + auxRight.size());
    res.insert(res.end(), auxLeft.begin(), auxLeft.end());
    res.insert(res.end(), auxRight.begin(), auxRight.end());

    return res;
}

std::vector<double> sh(const std::vector<double>& vec, int l, int r) {
    std::vector<double> aux, res;
    res = vec;
    int m = (l + r) / 2;

    if (res.size() % 2 != 0) {
        int lValue = res.at(res.size() - 1);

        res.erase(res.end() - 1);

        for (int j = l; j < (r - 1) / 2; j++) {
            aux.push_back(res.at(j));
            aux.push_back(res.at(j + m));
        }

        aux.push_back(lValue);
    } else {
        for (int j = l; j < r / 2; j++) {
            aux.push_back(res.at(j));
            aux.push_back(res.at(j + m));
        }
    }

    res.erase(res.begin(), res.end());
    res.shrink_to_fit();

    for (int i = l; i < r; i++) {
        res.push_back(aux[i]);
    }
    res = loc_sort(res, res.size());

    return res;
}

std::vector<double> BatcherSortTBB(const std::vector<double>& vec) {
    std::vector<double> res(vec);
    std::vector<double> auxLeft, auxRight;
    int size = vec.size();
    int half = size / 2;

    if (vec.size() == 2) loc_sort(vec, vec.size());
    if (vec.size() == 1) return vec;

    res = unsTBB(res, 0, size);

    tbb::task_scheduler_init init;

    tbb::parallel_invoke([&auxLeft, &res]
    { auxLeft = loc_sort(res, res.size()); },
        [&auxRight, &res] { auxRight = loc_sort(res, res.size()); });


    tbb::parallel_for(tbb::blocked_range<int>(0, half),
        [&res, &auxLeft](const tbb::blocked_range<int>& r) {
        int begin = r.begin();
        int end = r.end();
        for (int i = begin; i < end; i++) {
            res[i] = auxLeft[i];
        }
    }, tbb::simple_partitioner());

    tbb::parallel_for(tbb::blocked_range<int>(half, size),
        [&res, &auxRight](const tbb::blocked_range<int>& t) {
        int begin = t.begin();
        int end = t.end();
        for (int i = begin; i < end; i++) {
            res[i] = auxRight[i];
        }
    }, tbb::simple_partitioner());


    res = sh(res, 0, size);

    return res;
}

std::vector<double> unsTBB(const std::vector<double>& vec, int l, int r) {
    std::vector<double> auxLeft, auxRight, res;
    res = vec;
    if (res.size() % 2 != 0) {
        int lValue = res.at(res.size() - 1), lPos = res.size() - 1;

        res.erase(res.end() - 1);
        for (int j = l; j < r - 1; j += 2) {
            auxLeft.push_back(res.at(j));
            auxRight.push_back(res.at(j + 1));
        }

        if (lPos % 2 == 0)
            auxLeft.push_back(lValue);
        else
            auxRight.push_back(lValue);
    } else {
        for (int j = l; j < r; j += 2) {
            auxLeft.push_back(res.at(j));
            auxRight.push_back(res.at(j + 1));
        }
    }

    res.erase(res.begin(), res.end());
    res.shrink_to_fit();

    res.reserve(auxLeft.size() + auxRight.size());
    res.insert(res.end(), auxLeft.begin(), auxLeft.end());
    res.insert(res.end(), auxRight.begin(), auxRight.end());

    return res;
}
