// Copyright 2021 Dmitry Ryazanov
#include <tbb/tbb.h>
#include <gtest/gtest.h>
#include <omp.h>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include "./shell_and_batcher_sort_tbb.h"

TEST(getRandomVector, no_exceptions) { ASSERT_NO_THROW(getRandomVector(100)); }

TEST(getRandomVector, wrong_size) { ASSERT_ANY_THROW(getRandomVector(-100)); }

TEST(getRandomVector, works) {
  auto vect = getRandomVector(100);
  ASSERT_EQ(static_cast<int>(vect.size()), 100);
}

TEST(Shell_and_Batcher_sort, no_exceptions) {
  auto vect = getRandomVector(100);
  ASSERT_NO_THROW(SortShell(vect));
}

TEST(Shell_and_Batcher_sort, works_size_two) {
  auto vect = std::vector<int>({2, 1});
  print(vect);
  auto res = SortShell(vect);
  print(res);
  ASSERT_EQ(res, std::vector<int>({1, 2}));
}

TEST(Shell_and_Batcher_sort, works_size_4) {
  std::vector<int> vect = getRandomVector(4);
  print(vect);
  auto res = BatcherSort(vect);
  print(res);
}

TEST(Shell_and_Batcher_sort, works_size_10_OMP) {
  std::vector<int> vect = getRandomVector(10);
  print(vect);
  auto res = BatcherSortTBB(vect);
  print(res);
}

TEST(Shell_and_Batcher_sort, works_random_size_100) {
  auto vect = getRandomVector(100);
  auto sorted = SortShell(vect);
  std::sort(vect.begin(), vect.end());
  ASSERT_EQ(vect, sorted);
}

TEST(Shell_and_Batcher_sort, test_time) {
  auto vect = getRandomVector(1000);
  tbb::tick_count start = tbb::tick_count::now();
  auto res = BatcherSort(vect);
  tbb::tick_count end = tbb::tick_count::now();
  std::cout << "Seq sort time: " << (end - start).seconds() << std::endl;

  tbb::tick_count start2 = tbb::tick_count::now();
  auto resOMP = BatcherSortTBB(vect);
  tbb::tick_count end2 = tbb::tick_count::now();
  std::cout << "OMP sort time: " << (end2 - start2).seconds() << std::endl;

  ASSERT_EQ(res, resOMP);
}

TEST(Shell_and_Batcher_sort, works_random_vec_Shell_and_Batcher) {
  auto vect = getRandomVector(100);
  auto res1 = SortShell(vect);
  auto res2 = BatcherSort(vect);
  std::sort(vect.begin(), vect.end());
  ASSERT_EQ(res1, res2);
}

TEST(Shell_and_Batcher_sort, works_random_vec_Shell_and_Batcher_OMP) {
  auto vect = getRandomVector(10);
  auto res2 = BatcherSortTBB(vect);
  std::sort(vect.begin(), vect.end());
  ASSERT_EQ(vect, res2);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
