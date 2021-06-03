// Copyright 2021 Schelyanskova Anastasiia
#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <vector>

#include "./radix_sort_omp.h"

TEST(getRandomVector, wrong_size) { ASSERT_ANY_THROW(getRandomVector(-100)); }

TEST(getRandomVector, no_exceptions) { ASSERT_NO_THROW(getRandomVector(100)); }

TEST(getRandomVector, works) {
  std::vector<int> vect = getRandomVector(100);
  ASSERT_EQ(static_cast<int>(vect.size()), 100);
}

TEST(get_max_power, no_exceptions) {
  std::vector<int> vect = {1, 75, 87, 9, 400, 5};
  int power = 3;
  int max_power = get_max_power(vect);
  ASSERT_EQ(power, max_power);
}

TEST(radix_sort_omp, correct_works) {
  std::vector<int> vect = {1, 7, 8, 9, 4, 5};
  std::vector<int> correctAnswer = {1, 4, 5, 7, 8, 9};
  std::vector<int> myAnswer = RadixSort(vect, 6);
  ASSERT_EQ(correctAnswer, myAnswer);
}

// TEST(simple_merge, simple_merge) {
//  std::vector<int> vect = {2, 9, 11, 7, 13, 4, 12, 8, 1, 15, 22};
//  std::vector<int> sorted = RadixSortParallel(vect, 1, 4);
//  std::vector<int> correctAnswer = {1, 2, 4, 7, 8, 9, 11, 12, 13, 15, 22};
//  ASSERT_EQ(sorted, correctAnswer);
// }

int main(int argc, char **argv) {
  // int size = 10000000;
  // std::vector<int> unsorted = getRandomVector(size);
  //  for (int i = 0; i < size; i++) {
  //  std::cout << unsorted[i] << " ";
  // }
  // std::cout << "\n";*/
  // double t1 = omp_get_wtime();
  // std::vector<int> sorted = RadixSort(unsorted, size);
  // double t2 = omp_get_wtime();
  // std::vector<int> sortedParallel = RadixSortParallel(unsorted, size, 4);
  // double t3 = omp_get_wtime();
  // std::cout << "Sequential time: " << t2 - t1 << "\n";
  // std::cout << "Parallel time: " << t3 - t2 << "\n";
  //  for (int i = 0; i < size; i++) {
  //  std::cout << sorted[i] << " ";
  // }
  // std::cout << "\n";
  // for (int i = 0; i < sortedParallel.size(); i++) {
  //  std::cout << sortedParallel[i] << " ";
  // }
  // std::cout << "\n";
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
