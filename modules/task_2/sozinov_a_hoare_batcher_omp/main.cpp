// Copyright 2021 Sozinov Alex

#include <gtest/gtest.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "./hoare_batcher_omp.h"

TEST(GetRandomVector, get_vector) {
  ASSERT_NO_THROW(GetRandomVector(100));
}

TEST(GetRandomVector, vector_have_correct_size) {
  std::vector<double> vect = GetRandomVector(100);
  ASSERT_EQ(static_cast<int>(vect.size()), 100);
}

TEST(Sort, can_sort_one_thread) {
  std::vector<double> vect = GetRandomVector(100);
  ASSERT_NO_THROW(ParSort(&vect, 1));
}

TEST(Sort, can_sort_four_thread) {
  std::vector<double> vect = GetRandomVector(100);
  ASSERT_NO_THROW(ParSort(&vect, 4));
}

TEST(Sort, correct_sort_10000) {
  std::vector<double> vect = GetRandomVector(100);
  std::vector<double> copy = vect;
  std::vector<double> par_vect = vect;
  std::sort(copy.begin(), copy.end());
  SeqSort(&vect);

  ParSort(&par_vect, 1);

  ASSERT_EQ(copy, vect);
  ASSERT_EQ(copy, par_vect);
}

TEST(Sort, correct_sort_900000) {
  std::vector<double> vect = GetRandomVector(900);
  std::vector<double> copy = vect;
  std::vector<double> par_vect = vect;
  std::sort(copy.begin(), copy.end());
  SeqSort(&vect);

  ParSort(&par_vect, 4);

  ASSERT_EQ(copy, vect);
  ASSERT_EQ(copy, par_vect);
}

TEST(Sort, correct_sort_100000) {
  std::vector<double> vect = GetRandomVector(100);
  std::vector<double> copy = vect;
  std::vector<double> par_vect = vect;
  std::sort(copy.begin(), copy.end());
  SeqSort(&vect);

  ParSort(&par_vect, 3);

  ASSERT_EQ(copy, vect);
  ASSERT_EQ(copy, par_vect);
}

TEST(Sort, correct_sort_3000000) {
  std::vector<double> vect = GetRandomVector(300);
  std::vector<double> copy = vect;
  std::vector<double> par_vect = vect;
  std::sort(copy.begin(), copy.end());
  SeqSort(&vect);

  ParSort(&par_vect, 4);

  ASSERT_EQ(copy, vect);
  ASSERT_EQ(copy, par_vect);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
