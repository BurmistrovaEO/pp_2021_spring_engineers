// Copyright 2021 Schelyanskova Anastasiia
#include <gtest/gtest.h>
#include <tbb/tbb.h>
#include <vector>
#include <algorithm>
#include "./bitwise_sort_simple_tbb.h"


TEST(RadixSort_tbb, not_exceptions) {
    auto vect = getRandomVector(10);
    ASSERT_NO_THROW(RadixSort_tbb(vect));
}

TEST(RadixSort_tbb, one_not_exceptions) {
    auto vect = getRandomVector(1);
    ASSERT_NO_THROW(RadixSort_tbb(vect));
}

TEST(RadixSort_tbb, two_not_exceptions) {
    auto vect = getRandomVector(2);
    ASSERT_NO_THROW(RadixSort_tbb(vect));
}

TEST(RadixSort_tbb, two) {
    auto vect = std::vector<int>({5, 1});
    auto sorted = RadixSort_tbb(vect);
    ASSERT_EQ(sorted, std::vector<int>({1, 5}));
}

TEST(RadixSort_tbb, three) {
    auto vect = std::vector<int>({6, 4, 1});
    auto sorted = RadixSort_tbb(vect);
    ASSERT_EQ(sorted, std::vector<int>({1, 4, 6}));
}

TEST(RadixSort_tbb, five) {
    auto vect = std::vector<int>({3, 2, 1, 4, 5});
    auto sorted = RadixSort_tbb(vect);
    ASSERT_EQ(sorted, std::vector<int>({1, 2, 3, 4, 5}));
}

TEST(RadixSort_tbb, zero) {
    auto vect = std::vector<int>();
    auto sorted = RadixSort_tbb(vect);
    ASSERT_EQ(sorted, std::vector<int>());
}

TEST(RadixSort_tbb, num_100) {
    auto vect = getRandomVector(100);
    auto sorted = RadixSort_tbb(vect);
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, sorted);
}

TEST(RadixSort_tbb, num_1000) {
    auto vect = getRandomVector(1000);
    auto sorted = RadixSort_tbb(vect);
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, sorted);
}

TEST(RadixSort_tbb, num_10000) {
    auto vect = getRandomVector(10000);
    auto sorted = RadixSort_tbb(vect);
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, sorted);
}

TEST(RadixSort_tbb, get_time) {
    auto vect = getRandomVector(1000000);

    tbb::tick_count start_1 = tbb::tick_count::now();
    auto p_sorted = RadixSort_tbb(vect);
    tbb::tick_count end_1 = tbb::tick_count::now();
    std::cout << "Parallel sort time: " << (end_1 - start_1).seconds() << std::endl;

    tbb::tick_count start_2 = tbb::tick_count::now();
    auto l_sorted = RadixSort(vect, vect.size());
    tbb::tick_count end_2 = tbb::tick_count::now();
    std::cout << "Linear sort time: " << (end_2 - start_2).seconds() << std::endl;

    ASSERT_EQ(p_sorted, l_sorted);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
