// Copyright 2021 Gavrilov Nikita
#define _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING

#include <gtest/gtest.h>
#include <omp.h>
#include <vector>
#include <iostream>
#include <ctime>
#include <random>
#include "../../modules/task_2/gavrilov_n_bms_omp/BatcherMergingSort.h"

const int dataCount = 1000;
const int maxSplitCount = 6;

bool IsSorted(const std::vector<int>& data) {
    if (data.size() == 0)
        return true;
    for (size_t i = 0; i < data.size() - 1; i++) {
        if (data[i] > data[i + 1])
            return false;
    }
    return true;
}
std::vector<int> GetRandomVector(size_t size) {
    std::mt19937 rnd;
    rnd.seed(static_cast<unsigned int>(time(0)));

    std::vector<int> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = rnd() % 20000 - 10000;
    }

    return result;
}

TEST(Parallel_Operations_MPI, No_Throw_SplitCount1toMax) {
    for (int i = 1; i <= maxSplitCount; i++) {
        std::vector<int> data = GetRandomVector(dataCount);
        ASSERT_NO_THROW(Sort(&data, i));
    }
}

TEST(Parallel_Operations_MPI, Sorting_Right_SplitCount1toMax) {
    for (int i = 1; i <= maxSplitCount; i++) {
        std::vector<int> data = GetRandomVector(dataCount);
        Sort(&data, i);
        ASSERT_TRUE(IsSorted(data));
    }
}

TEST(Parallel_Operations_MPI, Sorting_Time_SplitCount1toMax) {
    for (int i = 1; i <= maxSplitCount; i++) {
        std::vector<int> data = GetRandomVector(dataCount);

        double begin = omp_get_wtime();
        Sort(&data, i);
        double end = omp_get_wtime();

        double t = end - begin;
        std::cout << "time of (splitCount = " << i << " count = " << dataCount << ") : " << t << std::endl;
    }
}

TEST(Parallel_Operations_MPI, Sort_By_Num_Place_Rang10) {
    std::vector<int> data = GetRandomVector(dataCount);
    ASSERT_NO_THROW(SortByNumPlace(&data, 10));
    ASSERT_TRUE(IsSorted(data));
}

TEST(Parallel_Operations_MPI, Sort_By_Num_Place_Rang2) {
    std::vector<int> data = GetRandomVector(dataCount);
    ASSERT_NO_THROW(SortByNumPlace(&data, 2));
    ASSERT_TRUE(IsSorted(data));
}

TEST(Parallel_Operations_MPI, Sorting_Right_Even_If_Bad_Size_SplitCount1toMax) {
    for (int j = 1; j <= maxSplitCount; j++) {
        for (int i = 1; i < j; i++) {
            std::vector<int> data = GetRandomVector(j * 1000 + j - i);
            ASSERT_NO_THROW(Sort(&data, j));
            ASSERT_TRUE(IsSorted(data));
        }
    }
}

TEST(Parallel_Operations_MPI, Sorting_Right_Even_If_Bad_Size2_SplitCount1toMax) {
    for (int j = 1; j <= maxSplitCount; j++) {
        for (int i = 1; i <= j; i++) {
            std::vector<int> data = GetRandomVector(j - 1);
            ASSERT_NO_THROW(Sort(&data, j));
            ASSERT_TRUE(IsSorted(data));
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
