// Copyright 2021 Bulychev Andrey
#include <omp.h>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./shell_sort_simple_omp.h"


TEST(Parallel_OMP, correct_work_parallel_sort) {
    std::vector<double> vec = getRandomVector(100);
    ASSERT_NO_THROW(Shell_sort_OMP(vec));
}

TEST(Parallel_OMP, work_parallel_sort_with_size_one) {
    std::vector<double> vec = getRandomVector(1);
    ASSERT_NO_THROW(Shell_sort_OMP(vec));
}

TEST(Parallel_OMP, compare_seq_and_parall) {
    std::vector<double> vec = getRandomVector(100);
    std::vector<double> seq = Shell_sort(vec);
    std::vector<double> omp = Shell_sort_OMP(vec);
    ASSERT_EQ(seq, omp);
}

TEST(Parallel_OMP, compare_seq_and_parall_with_bug_size) {
    std::vector<double> vec = getRandomVector(100);
    std::vector<double> seq = Shell_sort(vec);
    std::vector<double> omp = Shell_sort_OMP(vec);
    ASSERT_EQ(seq, omp);
}

TEST(Parallel_OMP, check_time) {
    std::vector<double> vec = getRandomVector(100);

    double start = omp_get_wtime();
    std::vector<double> omp = Shell_sort_OMP(vec);
    double end = omp_get_wtime();
    std::cout << "Parallel sort time: " << end - start << std::endl;

    start = omp_get_wtime();
    std::vector<double> seq = Shell_sort(vec);
    end = omp_get_wtime();
    std::cout << "Sequent sort time: " << end - start << std::endl;

    ASSERT_EQ(seq, omp);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
