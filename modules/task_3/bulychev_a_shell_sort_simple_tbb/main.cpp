// Copyright 2021 Bulychev Andrey
#include <tbb/tbb.h>
#include <omp.h>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./shell_sort_simple_tbb.h"


TEST(Parallel_OMP, correct_work_parallel_sort) {
  std::vector<double> vec = getRandomVector(100);
  ASSERT_NO_THROW(Shell_Sort_tbb(vec));
}

TEST(Parallel_TBB, work_parallel_sort_with_size_one) {
  std::vector<double> vec = getRandomVector(1);
  ASSERT_NO_THROW(Shell_Sort_tbb(vec));
}

TEST(Parallel_TBB, compare_seq_and_parall) {
  std::vector<double> vec = getRandomVector(100);
  std::vector<double> seq = Shell_sort(vec);
  std::vector<double> omp = Shell_Sort_tbb(vec);
  ASSERT_EQ(seq, omp);
}

TEST(Parallel_TBB, compare_seq_and_parall_with_bug_size) {
  std::vector<double> vec = getRandomVector(100);
  std::vector<double> seq = Shell_sort(vec);
  std::vector<double> omp = Shell_Sort_tbb(vec);
  ASSERT_EQ(seq, omp);
}

TEST(Parallel_TBB, check_time) {
  std::vector<double> vec = getRandomVector(1000);

  tbb::tick_count start = tbb::tick_count::now();
  std::vector<double> omp = Shell_Sort_tbb(vec);
  tbb::tick_count end = tbb::tick_count::now();
  std::cout << "Parallel sort time: " << (end - start).seconds() << std::endl;

  start = tbb::tick_count::now();
  std::vector<double> seq = Shell_sort(vec);
  end = tbb::tick_count::now();
  std::cout << "Sequent sort time: " << (end - start).seconds() << std::endl;

  ASSERT_EQ(seq, omp);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
