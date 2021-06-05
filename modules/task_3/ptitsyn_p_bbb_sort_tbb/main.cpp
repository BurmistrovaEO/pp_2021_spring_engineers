// Copyright 2018 Nesterov Alexander
#include <gtest/gtest.h>
#include <omp.h>
#include <vector>
#include "./bbb_sort.h"

TEST(RadixSortBatcherMerge, radix_one_element) {
    auto vect = random_double_array(1u);

    auto res1 = BatcherSort(vect);
    std::sort(vect.begin(), vect.end());

    ASSERT_EQ(vect, res1);
}

TEST(RadixSortBatcherMerge, radix_few_elements) {
    auto vect = random_double_array(5u);

    auto res1 = BatcherSort(vect);
    std::sort(vect.begin(), vect.end());

    ASSERT_EQ(vect, res1);
}

TEST(RadixSortBatcherMerge1, radix_lots_of_elements2) {
    auto vect = random_double_array(1000u);

    auto res1 = BatcherSort(vect);
    std::sort(vect.begin(), vect.end());

    ASSERT_EQ(vect, res1);
}

TEST(RadixSortBatcherMerge2, radix_lots_of_elements) {
    auto vect = random_double_array(10000u);

    auto res1 = BatcherSort(vect);
    std::sort(vect.begin(), vect.end());

    ASSERT_EQ(vect, res1);
}



TEST(Radix_Sort_W_Batcher, Test_Batcher) {
    auto vect = random_double_array(100u);
    double start_time = omp_get_wtime();
    auto res2 = BatcherSortTBB(vect);
    double end_time = omp_get_wtime();
    std::cout << "Sort: " << end_time - start_time << std::endl;
    start_time = omp_get_wtime();
    auto res1 = BatcherSort(vect);
    end_time = omp_get_wtime();
    std::cout << "Sort: " << end_time - start_time << std::endl;
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, res2);
    ASSERT_EQ(vect, res1);
}

TEST(Bitwise_and_Batcher_sort, works_random_vec_bitwise_Batcher_TBB) {
    auto vect = random_double_array(1000u);
    double start_time = omp_get_wtime();
    auto res2 = BatcherSortTBB(vect);
    double end_time = omp_get_wtime();
    std::cout << "Sort: " << end_time - start_time << std::endl;
    start_time = omp_get_wtime();
    auto res1 = BatcherSort(vect);
    end_time = omp_get_wtime();
    std::cout << "Sort: " << end_time - start_time << std::endl;
    std::sort(vect.begin(), vect.end());
    ASSERT_EQ(vect, res2);
    ASSERT_EQ(vect, res1);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
