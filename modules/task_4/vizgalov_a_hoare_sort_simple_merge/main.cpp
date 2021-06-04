// Copyright 2021 Vizgalov Anton

#include <gtest/gtest.h>
#include <algorithm>

#include "../../../3rdparty/unapproved/unapproved.h"
#include "./hoare_sort_simple_merge.h"

TEST(Hoare_Sort_Simple_Merge_STD, Hoare_Sort_Empty_Vector) {
    std::vector<int> vec = getRandomVector(0);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoareSort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Hoare_Sort_Only_Positive_Numbers) {
    std::vector<int> vec = getRandomVector(100, 1);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoareSort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Hoare_Sort_All_Numbers_Small_Size) {
    std::vector<int> vec = getRandomVector(100);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoareSort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Hoare_Sort_All_Numbers_Big_Size) {
    std::vector<int> vec = getRandomVector(10000);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoareSort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Hoare_Sort_Already_Sorted) {
    std::vector<int> vec = getRandomVector(100, 2);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoareSort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Hoare_Sort_Worst_Case_Unsorted) {
    std::vector<int> vec = getRandomVector(100, 3);
    int lastIdx = static_cast<int>(vec.size() - 1);

    ASSERT_NO_THROW(hoareSort(&vec, 0, lastIdx));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Merge_Left_Empty) {
    std::vector<int> vec = getRandomVector(100);
    int lastIdx = static_cast<int>(vec.size() - 1);
    hoareSort(&vec, 0, lastIdx);

    ASSERT_NO_THROW(merge(&vec, 0, 0, vec.size() - 1));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Merge_Right_Empty) {
    std::vector<int> vec = getRandomVector(100);
    int lastIdx = static_cast<int>(vec.size() - 1);
    hoareSort(&vec, 0, lastIdx);

    ASSERT_NO_THROW(merge(&vec, 0, vec.size() - 1, vec.size() - 1));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Merge_Same_Size) {
    std::vector<int> vec = getRandomVector(100);
    int lastIdx = static_cast<int>(vec.size() - 1);
    int midIdx = vec.size() / 2;
    hoareSort(&vec, 0, midIdx);
    hoareSort(&vec, midIdx + 1, lastIdx);

    ASSERT_NO_THROW(merge(&vec, 0, midIdx, vec.size() - 1));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Merge_Different_Size) {
    std::vector<int> vec = getRandomVector(100);
    int lastIdx = static_cast<int>(vec.size() - 1);
    int divisionIdx = vec.size() / 8;
    hoareSort(&vec, 0, divisionIdx);
    hoareSort(&vec, divisionIdx + 1, lastIdx);

    ASSERT_NO_THROW(merge(&vec, 0, divisionIdx, vec.size() - 1));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Merge_Large_Vectors) {
    std::vector<int> vec = getRandomVector(10000);
    int lastIdx = static_cast<int>(vec.size() - 1);
    int midIdx = vec.size() / 2;
    hoareSort(&vec, 0, midIdx);
    hoareSort(&vec, midIdx + 1, lastIdx);

    ASSERT_NO_THROW(merge(&vec, 0, midIdx, vec.size() - 1));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Run_Hoare_Sort_Even_Size) {
    std::vector<int> vec = getRandomVector(1000000);

    ASSERT_NO_THROW(runHoareSort(&vec));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Run_Hoare_Sort_Odd_Size) {
    std::vector<int> vec = getRandomVector(1111111);

    ASSERT_NO_THROW(runHoareSort(&vec));

    EXPECT_TRUE(std::is_sorted(vec.begin(), vec.end()));
}

TEST(Hoare_Sort_Simple_Merge_STD, Compare_Time) {
  std::vector<int> vec1 = getRandomVector(1000000);

  auto t0 = std::chrono::steady_clock::now();
  ASSERT_NO_THROW(runHoareSort(&vec1));
  auto t1 = std::chrono::steady_clock::now();
  std::cout << "Time (parallel): " << ((t1 - t0).count() / 1e6) << "ms" << std::endl;

  std::vector<int> vec2 = getRandomVector(1000000);

  auto t2 = std::chrono::steady_clock::now();
  ASSERT_NO_THROW(hoareSort(&vec2, 0, vec2.size() - 1));
  auto t3 = std::chrono::steady_clock::now();
  std::cout << "Time (non-parallel): " << ((t3 - t2).count() / 1e6) << "ms" << std::endl;

  EXPECT_TRUE(std::is_sorted(vec1.begin(), vec1.end()));
  EXPECT_TRUE(std::is_sorted(vec2.begin(), vec2.end()));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
