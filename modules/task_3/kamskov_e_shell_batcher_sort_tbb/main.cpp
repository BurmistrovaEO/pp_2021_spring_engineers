// Copyright 2021 Kamskov Eugene
#include <gtest/gtest.h>
#include <omp.h>
#include <tbb/tbb.h>
#include "./shell_batcher.h"

TEST(test_Shell, shell_Rand) {
    std::vector<int> test = getRandomVector(20);
    std::vector<int> shell = shellSort(test, 20);
    std::sort(test.begin(), test.end());
    ASSERT_EQ(shell, test);
}

TEST(test_Shell_Betch_sort_Tbb, test_on_1000_elements_4_threads) {
    std::vector<int> res = getRandomVector(1001);
    std::vector<int> betch = batchShellSortTbb(res, 4, 1001);
    std::sort(res.begin(), res.end());
    ASSERT_EQ(betch, res);
}
TEST(test_Batcher_Shell_Dif_Size_Tbb, test_time_4_threads) {
    std::vector<int> vect = getRandomVector(1000001);
    tbb::tick_count start = tbb::tick_count::now();
    std::vector<int> res = batchShellSort(vect, 1000001);
    tbb::tick_count end = tbb::tick_count::now();
    std::cout << "Seq time: " << (end - start).seconds() << std::endl;

    tbb::tick_count start2 = tbb::tick_count::now();
    std::vector<int> resOMP = batchShellSortTbb(vect, 4, 1000001);
    tbb::tick_count end2 = tbb::tick_count::now();
    std::cout << "Tbb time: " << (end2 - start2).seconds() << std::endl;

    ASSERT_EQ(res, resOMP);
}

TEST(test_Batcher_Manual_Dif_Size, batcher_Manual) {
    std::vector<int> vec1 = { 1, 10, 15, 20, 30, 11, 2 };
    std::vector<int> vec2 = { 0, 3, 5, 7, 40 };
    std::sort(vec1.begin(), vec1.end());
    std::sort(vec2.begin(), vec2.end());
    std::vector<int>stl_res(12);
    std::merge(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), stl_res.begin());
    std::vector<int>even = evenBatch(vec1, vec2);
    std::vector<int>odd = oddBatch(vec1, vec2);
    std::vector<int>my_res = mergeBatch(even, odd);
    ASSERT_EQ(stl_res, my_res);
}

TEST(test_Batcher_Shell_Tbb, batcher_shell_Tbb) {
    std::vector<int> vect = getRandomVector(1000);
    std::vector<int> res2 = batchShellSortTbb(vect, 4, 1000);
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, res2);
}

TEST(test_Shell_Batch_Manual_Tbb, shell_Batch_Manual_Omp_4_threads) {
    std::vector<int> stl_res = { 1, 10, 15, 20, 30, 11, 2, 32, 11, 432, 222, 6453, 76, 3 };
    std::vector<int> my_res = batchShellSortTbb(stl_res, 4, 14);
    std::sort(stl_res.begin(), stl_res.end());
    ASSERT_EQ(my_res, stl_res);
}

TEST(test_Shell_Batch_Rand, shell_Batch_Rand) {
    std::vector<int> res = getRandomVector(100);
    std::vector<int> my_res = batchShellSortTbb(res, 4, 100);
    std::sort(res.begin(), res.end());
    ASSERT_EQ(my_res, res);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
