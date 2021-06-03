// Copyright 2018 Nesterov Alexander
#include <gtest/gtest.h>
#include <vector>
#include "./lin_hist_stret.h"

TEST(Parallel_Operations_TBB, Test_random_matrix_with_negative_size) {
    std::vector<int> empty;  // empty vector
    int size = -5;
    ASSERT_EQ(getRandomMatrix(size, size), empty);
}

TEST(Parallel_Operations_TBB, Test_negative_size_stretching) {
    std::vector<int> empty;  // empty vector
    std::vector<int> test;
    int size = -5;
    ASSERT_EQ(histogramStretching(test, size, size), empty);
}

TEST(Parallel_Operations_TBB, Test_values_range_of_random_matrix) {
    int size = 4;
    int rightRange = 1;
    std::vector<int> test = getRandomMatrix(size, size);
    for (int i = 0; i < size * size; i++) {
        if (test[i] < 0 || test[i]>255) {
            rightRange = 0;
        }
    }
    ASSERT_EQ(rightRange, 1);
}

TEST(Parallel_Operations_TBB, Test_stretching_of_manual_matrix) {
    int size = 3;
    std::vector<int> test{ 100, 245, 100,
                           45, 245, 45,
                           100, 245, 100 };
    std::vector<int> result{ 55, 200, 55,
                             0, 200, 0,
                             55, 200, 55 };
    ASSERT_EQ(histogramStretching(test, size, size), result);
}

TEST(Parallel_Operations_TBB, Test_stretching_of_random_matrix) {
    int size = 1000;
    std::vector<int> test = getRandomMatrix(size, size);
    std::vector<int> empty;  // empty vector
    test = histogramStretching(test, size, size);
    ASSERT_NE(test, empty);
}

TEST(Parallel_Operations_TBB, Test_negative_size_tbb_stretching) {
    std::vector<int> empty;  // empty vector
    std::vector<int> test;
    int size = -5;
    ASSERT_EQ(histogramStretchingTbb(test, size, size), empty);
}

TEST(Parallel_Operations_TBB, Comparison_seq_with_tbb) {
    int size = 5000;
    std::vector<int> mx = getRandomMatrix(size, size);
    std::vector<int> seq;
    std::vector<int> tbb;
    seq = histogramStretching(mx, size, size);
    tbb = histogramStretchingTbb(mx, size, size);
    ASSERT_EQ(seq, tbb);
}

TEST(Parallel_Operations_TBB, Comparison_hists) {
    int size = 5000;
    std::vector<int> mx = getRandomMatrix(size, size);
    std::vector<int> seq;
    std::vector<int> tbb;
    seq = histogramStretching(mx, size, size);
    tbb = histogramStretchingTbb(mx, size, size);
    std::vector<int> seq_hist(256, 0);
    std::vector<int> tbb_hist(256, 0);
    for (int i = 1; i < size*size; i++) {
        seq_hist[seq[i]]++;
        tbb_hist[tbb[i]]++;
    }
    ASSERT_EQ(seq_hist, tbb_hist);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
