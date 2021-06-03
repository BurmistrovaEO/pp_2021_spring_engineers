// Copyright 2021 Kamskov Eugene

#include <gtest/gtest.h>
#include "./shell_batcher.h"

TEST(test_Shell, shell_Rand) {
    std::vector<int> test = getRandomVector(20);
    std::vector<int> shell = shellSort(test, 20);
    std::sort(test.begin(), test.end());
    ASSERT_EQ(shell, test);
}

TEST(test_Shell_Betch_sort_omp, test_on_1000_elements_4_threads) {
    std::vector<int> res = getRandomVector(1001);
    std::vector<int> betch = batchShellSortOmp(res, 4, 1001);
    std::sort(res.begin(), res.end());
    ASSERT_EQ(betch, res);
}
TEST(test_Batcher_Shell_Dif_Size_Omp, test_time_10_threads) {
    std::vector<int> vect = getRandomVector(1000001);
    double start = omp_get_wtime();
    std::vector<int> res = batchShellSort(vect, 1000001);
    double end = omp_get_wtime();
    std::cout << "Seq sort time: " << end - start << std::endl;

    start = omp_get_wtime();
    std::vector<int> resOMP = batchShellSortOmp(vect, 10, 1000001);
    end = omp_get_wtime();
    std::cout << "OMP sort time: " << end - start << std::endl;

    ASSERT_EQ(res, resOMP);
}

TEST(test_Batcher_Manual_Dif_Size, batcher_Manual) {
    std::vector<int> vec1 = {1, 10, 15, 20, 30, 11, 2};
    std::vector<int> vec2 = {0, 3, 5, 7, 40};
    std::sort(vec1.begin(), vec1.end());
    std::sort(vec2.begin(), vec2.end());
    std::vector<int>stl_res(12);
    std::merge(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), stl_res.begin());
    std::vector<int>even = evenBatch(vec1, vec2);
    std::vector<int>odd = oddBatch(vec1, vec2);
    std::vector<int>my_res = mergeBatch(even, odd);
    ASSERT_EQ(stl_res, my_res);
}

TEST(test_Batcher_Shell_Omp, batcher_shell_Omp) {
    std::vector<int> vect = getRandomVector(1000);
    std::vector<int> res2 = batchShellSortOmp(vect, 4, 1000);
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, res2);
}

TEST(test_Shell_Batch_Manual_Omp, shell_Batch_Manual_Omp_4_threads) {
    std::vector<int> stl_res = {1, 10, 15, 20, 30, 11, 2, 32, 11, 432, 222, 6453, 76, 3};
    std::vector<int> my_res = batchShellSortOmp(stl_res, 4, 14);
    std::sort(stl_res.begin(), stl_res.end());
    ASSERT_EQ(my_res, stl_res);
}

TEST(test_Shell_Batch_Rand, shell_Batch_Rand) {
    std::vector<int> res = getRandomVector(100);
    std::vector<int> my_res = batchShellSortOmp(res, 4, 100);
    std::sort(res.begin(), res.end());
    ASSERT_EQ(my_res, res);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
