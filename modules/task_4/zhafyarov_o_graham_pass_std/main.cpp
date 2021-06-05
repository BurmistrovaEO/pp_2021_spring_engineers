// Copyright 2021 Zhafyarov Oleg
#include <gtest/gtest.h>
#include <vector>
#include "./graham_pass.h"

TEST(STD, x100_Points_Check) {
  int size = 100;
  std::vector<point> vec = RandomVector(size);
  std::vector<point> result_std;
  std::vector<point> result_seq;
  result_std = GrahamPassStd(vec);
  result_seq = GrahamPassSeq(vec);
  ASSERT_EQ(CompareVectors(result_std, result_seq), true);
}

TEST(STD, x200_Points_Check) {
  int size = 200;
  std::vector<point> vec = RandomVector(size);
  std::vector<point> result_std;
  std::vector<point> result_seq;
  result_std = GrahamPassStd(vec);
  result_seq = GrahamPassSeq(vec);
  ASSERT_EQ(CompareVectors(result_std, result_seq), true);
}

TEST(STD, Count_Check) {
  int count_tmp;
  std::vector<point> vec = { {1, 1}, {1, 4},
                                           {4, 1}, {4, 4}, {2, 2}};
  std::vector<point> result;
  count_tmp = 4;
  result = GrahamPassStd(vec);
  ASSERT_EQ(static_cast<int>(result.size()), count_tmp);
}

TEST(STD, x300_Points_Check) {
  int size = 300;
  std::vector<point> vec = RandomVector(size);
  std::vector<point> result_std;
  std::vector<point> result_seq;
  result_std = GrahamPassStd(vec);
  result_seq = GrahamPassSeq(vec);
  ASSERT_EQ(CompareVectors(result_std, result_seq), true);
}

TEST(STD, x500_Points_Check) {
  int size = 500;
  std::vector<point> vec = RandomVector(size);
  std::vector<point> result_std;
  std::vector<point> result_seq;
  result_std = GrahamPassStd(vec);
  result_seq = GrahamPassSeq(vec);
  ASSERT_EQ(CompareVectors(result_std, result_seq), true);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
