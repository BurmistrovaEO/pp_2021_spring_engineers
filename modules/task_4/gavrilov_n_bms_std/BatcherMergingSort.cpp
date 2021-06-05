// Copyright 2021 Gavrilov Nikita
#include <utility>
#include <algorithm>
#include <vector>
#include <cmath>
#include "../../3rdparty/unapproved/unapproved.h"
#include "../../modules/task_4/gavrilov_n_bms_std/BatcherMergingSort.h"

#define MY_INT_MAX 2147483647

int8_t GetDigit(int value, uint8_t digidNum, uint8_t rang) {
    if (rang <= 1)
        throw "rang must be more than 1";

    int powRang = static_cast<int>(pow(rang, digidNum));
    int powRangNext = static_cast<int>(pow(rang, digidNum + 1.0));

    return value % powRangNext / powRang;
}
uint8_t GetNumOfDigits(int value, uint8_t rang) {
    if (rang <= 1)
        throw "rang must be more than 1";

    int length = 1;
    while (value /= rang)
        length++;

    return length;
}

std::vector<int> ConcatVectors(const std::vector<std::vector<int>>& vectors) {
    std::vector<int> result;
    for (size_t i = 0; i < vectors.size(); i++) {
        for (size_t j = 0; j < vectors[i].size(); j++) {
            result.push_back(vectors[i][j]);
        }
    }
    return result;
}
void SortByDigid(std::vector<int>* data, uint8_t digitNum, uint8_t rang) {
    if (rang <= 1)
        throw "rang must be more than 1";
    std::vector<std::vector<int>> vectorsByDigits(rang * 2 - 1);

    for (size_t i = 0; i < data->size(); i++) {
        int8_t digit = GetDigit((*data)[i], digitNum, rang);
        vectorsByDigits[digit + rang - 1].push_back((*data)[i]);
    }

    *data = ConcatVectors(vectorsByDigits);
}
void SortByNumPlace(std::vector<int>* data, uint8_t rang) {
    if (rang <= 1)
        throw "rang must be more than 1";

    uint8_t maxNumOfDigits = 0;

    for (size_t i = 0; i < data->size(); i++) {
        uint8_t numOfDigits = GetNumOfDigits((*data)[i], rang);
        if (numOfDigits > maxNumOfDigits)
            maxNumOfDigits = numOfDigits;
    }

    for (uint8_t i = 0; i < maxNumOfDigits; i++) {
        SortByDigid(data, i, rang);
    }
}

void SplitEvenOdd(const std::vector<int>& from, std::vector<int>* destEven, std::vector<int>* destOdd) {
    for (size_t i = 0; i < from.size() / 2; i++) {
        destEven->push_back(from[i * 2]);
        destOdd->push_back(from[i * 2 + 1]);
    }
    if (from.size() % 2 == 1) {
        destEven->push_back(from.back());
    }
}
void BatcherMerge(std::vector<int> procsLeft, std::vector<int> procsRight, std::vector<std::pair<int, int>>* comps) {
    size_t procsCount = procsLeft.size() + procsRight.size();
    if (procsCount <= 1) {
        return;
    } else if (procsCount == 2) {
        comps->push_back(std::pair<int, int>(procsLeft[0], procsRight[0]));
        return;
    }

    std::vector<int> procsLeftOdd, procsLeftEven;
    std::vector<int> procsRightOdd, procsRightEven;

    SplitEvenOdd(procsLeft, &procsLeftEven, &procsLeftOdd);
    SplitEvenOdd(procsRight, &procsRightEven, &procsRightOdd);

    BatcherMerge(procsLeftOdd, procsRightOdd, comps);
    BatcherMerge(procsLeftEven, procsRightEven, comps);

    std::vector<std::vector<int>> vecs{ procsLeft, procsRight };
    std::vector<int> procsResult = ConcatVectors(vecs);

    for (size_t i = 1; i + 1 < procsResult.size(); i += 2) {
        comps->push_back(std::pair<int, int>(procsResult[i], procsResult[i + 1]));
    }
}
void BatcherSplitNMerge(std::vector<int> procs, std::vector<std::pair<int, int>>* comps) {
    if (procs.size() <= 1) {
        return;
    }

    std::vector<int> procs_up(procs.begin(), procs.begin() + procs.size() / 2);
    std::vector<int> procs_down(procs.begin() + procs.size() / 2, procs.end());

    BatcherSplitNMerge(procs_up, comps);
    BatcherSplitNMerge(procs_down, comps);
    BatcherMerge(procs_up, procs_down, comps);
}
std::vector<std::pair<int, int>> Batcher(int procCount) {
    std::vector<int> procs(procCount);
    for (int i = 0; i < procCount; i++) {
        procs[i] = i;
    }

    std::vector<std::pair<int, int>> comps;
    BatcherSplitNMerge(procs, &comps);

    return comps;
}

void Sort(std::vector<int>* data, int splitCount) {
    if (data->size() == 0)
        return;

    size_t realSize = data->size();

    int size = splitCount;

    int maxSize = size < static_cast<int>(data->size()) ? size : static_cast<int>(data->size());
    std::vector<std::pair<int, int>> comps = Batcher(maxSize);
    while (data->size() % size != 0) {
        data->push_back(MY_INT_MAX);
    }
    int countPerProc = static_cast<int>(data->size()) / size;

    std::vector<std::thread> threads;
    threads.reserve(splitCount);

    for (int i = 0; i < splitCount; i++) {
        threads.push_back(std::thread([i, &data, countPerProc] {
            auto b = data->begin() + i * countPerProc;
            auto e = b + countPerProc;
            std::vector<int> localData = std::vector<int>(b, e);
            SortByNumPlace(&localData);
            std::copy(localData.begin(), localData.end(), data->begin() + i * countPerProc);
        }));
    }

    for (size_t i = 0; i < threads.size(); i++) {
        threads[i].join();
    }

    if (size > 1) {
        for (size_t i = 0; i < comps.size(); i++) {
            int beginIndexOfFirst = comps[i].first * countPerProc;
            int beginIndexOfSecond = comps[i].second * countPerProc;


            std::vector<int> tmp;
            tmp.reserve(countPerProc * 2);
            int j = 0, k = 0;
            for (int l = 0; l < countPerProc * 2; l++) {
                int currentIndexOfSecond = beginIndexOfSecond + j;
                int currentIndexOfFirst = beginIndexOfFirst + k;
                int index;
                if (j < countPerProc &&
                    (k >= countPerProc || data->at(currentIndexOfSecond) < data->at(currentIndexOfFirst))) {
                    index = currentIndexOfSecond;
                    j++;
                } else {
                    index = currentIndexOfFirst;
                    k++;
                }
                tmp.push_back(data->at(index));
            }

            std::copy(tmp.begin(), tmp.begin() + countPerProc, data->begin() + beginIndexOfFirst);
            std::copy(tmp.begin() + countPerProc, tmp.end(), data->begin() + beginIndexOfSecond);
        }
    }

    data->erase(data->begin() + realSize, data->end());
}
