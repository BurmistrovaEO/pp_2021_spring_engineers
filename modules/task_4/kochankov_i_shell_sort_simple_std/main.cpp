// Copyright 2021 Kochankov Ilya
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>
#include "../../../3rdparty/unapproved/unapproved.h"
#include "./shell_sort_simple_std.h"

TEST(shell_sort_std, no_exceptions) {
  auto vect = getRandomVector(100);
  ASSERT_NO_THROW(shell_sort_std(vect));
}

TEST(shell_sort_std, size_one_no_exceptions) {
  auto vect = getRandomVector(1);
  ASSERT_NO_THROW(shell_sort_std(vect));
}

TEST(shell_sort_std, size_two_no_exceptions) {
  auto vect = getRandomVector(2);
  ASSERT_NO_THROW(shell_sort_std(vect));
}

TEST(shell_sort_std, works_size_two) {
  auto vect = std::vector<double>({2, 1});
  auto sorted = shell_sort_std(vect);
  ASSERT_EQ(sorted, std::vector<double>({1, 2}));
}

TEST(shell_sort_std, works_size_three) {
  auto vect = std::vector<double>({3, 2, 1});
  auto sorted = shell_sort_std(vect);
  ASSERT_EQ(sorted, std::vector<double>({1, 2, 3}));
}

TEST(shell_sort_std, works_size_five) {
  auto vect = std::vector<double>({3, 2, 1, 6, 7});
  auto sorted = shell_sort_std(vect);
  ASSERT_EQ(sorted, std::vector<double>({1, 2, 3, 6, 7}));
}

TEST(shell_sort_std, works_size_zero) {
  auto vect = std::vector<double>();
  auto sorted = shell_sort_std(vect);
  ASSERT_EQ(sorted, std::vector<double>());
}

TEST(shell_sort_std, works_random_size_100) {
  auto vect = getRandomVector(100);
  auto sorted = shell_sort_std(vect);
  std::sort(vect.begin(), vect.end());
  ASSERT_EQ(vect, sorted);
}

TEST(shell_sort_std, works_random_size_1000) {
  auto vect = getRandomVector(1000);
  auto sorted = shell_sort_std(vect);
  std::sort(vect.begin(), vect.end());
  ASSERT_EQ(vect, sorted);
}

TEST(shell_sort_std, works_random_size_10000) {
  auto vect = getRandomVector(10000);
  auto sorted = shell_sort_std(vect);
  std::sort(vect.begin(), vect.end());
  ASSERT_EQ(vect, sorted);
}

TEST(shell_sort_std, get_time) {
  auto vect = getRandomVector(1000000);

  auto t0 = std::chrono::steady_clock::now();
  auto p_sorted = shell_sort_std(vect);
  auto t1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_ms = t1 - t0;
  std::cout << "Parallel sort time: " << elapsed_ms.count() << std::endl;

  t0 = std::chrono::steady_clock::now();
  auto l_sorted = shell_sort(vect);
  t1 = std::chrono::steady_clock::now();
  elapsed_ms = t1 - t0;
  std::cout << "Linear sort time: " << elapsed_ms.count() << std::endl;

  ASSERT_EQ(p_sorted, l_sorted);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
