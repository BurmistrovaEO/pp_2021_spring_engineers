// Copyright 2021 Solomakhin Sergey

#include <tbb/tbb.h>
#include <algorithm>
#include "gtest/gtest.h"

#include "./hoar_sort_simple_merge.h"

TEST(Hoare_Sort_Simple_Merge_TBB, Hoare_Sort_Empty_Vector) {
    std::vector<int> vec = random_gen(0);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoar_sort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Hoare_Sort_Only_Positive_Numbers) {
    std::vector<int> vec = random_gen(100, 1);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoar_sort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Hoare_Sort_All_Numbers_Small_Size) {
    std::vector<int> vec = random_gen(100);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoar_sort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Hoare_Sort_All_Numbers_Big_Size) {
    std::vector<int> vec = random_gen(10000);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoar_sort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Hoare_Sort_Already_Sorted) {
    std::vector<int> vec = random_gen(100, 2);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoar_sort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Hoare_Sort_Worst_Case_Unsorted) {
    std::vector<int> vec = random_gen(100, 3);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoar_sort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Merge_Left_Empty) {
    std::vector<int> vec1 = random_gen(100);
    std::vector<int> vec2 = random_gen(0);
    hoar_sort(&vec1, 0, vec1.size() - 1);

    std::vector<int> vec = merge(vec1, vec2);

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Merge_Right_Empty) {
    std::vector<int> vec1 = random_gen(0);
    std::vector<int> vec2 = random_gen(100);
    hoar_sort(&vec2, 0, vec2.size() - 1);

    std::vector<int> vec = merge(vec1, vec2);

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Merge_Same_Size) {
    std::vector<int> vec1 = random_gen(100);
    std::vector<int> vec2 = random_gen(100);
    hoar_sort(&vec1, 0, vec1.size() - 1);
    hoar_sort(&vec2, 0, vec2.size() - 1);

    std::vector<int> vec = merge(vec1, vec2);

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Merge_Different_Size) {
    std::vector<int> vec1 = random_gen(75);
    std::vector<int> vec2 = random_gen(25);
    hoar_sort(&vec1, 0, vec1.size() - 1);
    hoar_sort(&vec2, 0, vec2.size() - 1);

    std::vector<int> vec = merge(vec1, vec2);

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Merge_Large_Vectors) {
    std::vector<int> vec1 = random_gen(5000);
    std::vector<int> vec2 = random_gen(5000);
    hoar_sort(&vec1, 0, vec1.size() - 1);
    hoar_sort(&vec2, 0, vec2.size() - 1);

    std::vector<int> vec = merge(vec1, vec2);

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Run_Hoare_Sort_Even_Size) {
    std::vector<int> vec = random_gen(100000);

    ASSERT_NO_THROW(hoar_sort_tbb(&vec));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Run_Hoare_Sort_Odd_Size) {
    std::vector<int> vec = random_gen(111111);

    ASSERT_NO_THROW(hoar_sort_tbb(&vec));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_TBB, Compare_Time) {
    std::vector<int> vec1 = random_gen(1000000);

    tbb::tick_count t0 = tbb::tick_count::now();
    ASSERT_NO_THROW(hoar_sort_tbb(&vec1, 4));
    tbb::tick_count t1 = tbb::tick_count::now();
    std::cout << "Time (parallel): " << (t1 - t0).seconds() << std::endl;

    std::vector<int> vec2 = random_gen(1000000);

    tbb::tick_count t3 = tbb::tick_count::now();
    ASSERT_NO_THROW(hoar_sort(&vec2, 0, vec2.size() - 1));
    tbb::tick_count t4 = tbb::tick_count::now();
    std::cout << "Time (non-parallel): " << (t4 - t3).seconds() << std::endl;
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
