// Copyright 2021 Kamskov Eugene

#include <gtest/gtest.h>
#include "./shell_batcher.h"

TEST(test_Rand, rand) {
    auto vect = getRandomVector(432);
    ASSERT_EQ(static_cast<int>(vect.size()), 432);
}

TEST(test_Rand_Wrong_Size, wrong_Size) {
    ASSERT_ANY_THROW(getRandomVector(-1));
}

TEST(test_Shell, shell_Rand) {
    std::vector<int> test = getRandomVector(20);
    std::vector<int> shell = shellSort(test, 20);
    std::sort(test.begin(), test.end());
    ASSERT_EQ(shell, test);
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

TEST(test_Batcher_Manual_Dif_Random_Size, batcher_Manual_Rand) {
    std::vector<int> vec1 = getRandomVector(20);
    std::vector<int> vec2 = getRandomVector(100);
    std::sort(vec1.begin(), vec1.end());
    std::sort(vec2.begin(), vec2.end());
    std::vector<int>stl_res(120);
    std::merge(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), stl_res.begin());
    std::vector<int>even = evenBatch(vec1, vec2);
    std::vector<int>odd = oddBatch(vec1, vec2);
    std::vector<int>my_res = mergeBatch(even, odd);
    ASSERT_EQ(stl_res, my_res);
}

TEST(test_Shell_Batch_Manual, shell_Batch_Manual) {
    std::vector<int> stl_res = {1, 10, 15, 20, 30, 11, 2, 32, 11, 432, 222, 6453, 76, 3};
    std::vector<int> my_res = batchShellSort(stl_res, 14);
    std::sort(stl_res.begin(), stl_res.end());
    ASSERT_EQ(my_res, stl_res);
}

TEST(test_Shell_Batch_Rand, shell_Batch_Rand) {
    std::vector<int> stl_res = getRandomVector(100);
    std::vector<int> my_res = batchShellSort(stl_res, 100);
    std::sort(stl_res.begin(), stl_res.end());
    ASSERT_EQ(my_res, stl_res);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
