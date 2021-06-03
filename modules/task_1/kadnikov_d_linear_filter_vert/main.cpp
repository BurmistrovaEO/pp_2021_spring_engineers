// Copyright 2021 Kadnikov Denis
#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "../../../modules/task_1/kadnikov_d_linear_filter_vert/linear_filter_vert.h"

TEST(Linear_filter_test_1, can_not_create_random_vector) {
  int r = -1;
  int c = -1;
  ASSERT_ANY_THROW(createRandomVector(r, c));
}

TEST(Linear_filter_test_2, gaussianCalculate_with_negative_sigma) {
  ASSERT_ANY_THROW(gaussianCalculate(5, -1));
}

TEST(Linear_filter_test_3, can_apply_1x1_vector) {
  int r = 1;
  int c = 1;

  std::vector<double> v(r * c);
  v[0] = 221;
  v = transposition(v, r, c);
  std::swap(c, r);
  std::vector<double> vector_calculation(r * c);
  vector_calculation[0] = 220;
  ASSERT_EQ(vector_calculation, gaussianFilterApply(v, r, c, 1, 5));
}

TEST(Linear_filter_test_4, can_apply_3x3_vector) {
  int r = 3;
  int c = 3;
  std::vector<double> vector(r * c);
  vector[0] = 1;
  vector[1] = 2;
  vector[2] = 3;
  vector[3] = 4;
  vector[4] = 5;
  vector[5] = 6;
  vector[6] = 7;
  vector[7] = 8;

  vector = transposition(vector, r, c);
  std::swap(c, r);
  std::vector<double> vector_calculation(r * c);

  vector_calculation[0] = 8;
  vector_calculation[1] = 7;
  vector_calculation[2] = 6;
  vector_calculation[3] = 5;
  vector_calculation[4] = 4;
  vector_calculation[5] = 3;
  vector_calculation[6] = 2;
  vector_calculation[7] = 1;

  ASSERT_NO_THROW(gaussianFilterApply(vector, r, c, 1, 5));
}

TEST(Linear_filter_test_5, gaussianCalculate_with_negative_radius) {
  ASSERT_ANY_THROW(gaussianCalculate(-5, 1));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
