// Copyright 2021 Ognev Denis

#include <gtest/gtest.h>

#include "../../task_3/ognev_d_markings_hulls_tbb/markings_hulls_tbb.h"

TEST(Sequential_Test, Marking_Test) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  false, true,  true,  true, true,
                      false, true,  false, false, false, false, true, true,
                      true,  true,  false, true,  false, false, true, false,
                      false, true,  true,  true,  false, true};
  std::vector<int> right{2, 0, 2, 0, 3, 2, 2, 2, 0, 3, 0, 0, 0, 0, 3,
                         4, 4, 4, 0, 3, 0, 0, 4, 0, 0, 4, 4, 4, 0, 6};
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (v[i * width + j] == true)
        std::cout << static_cast<char>(219);
      else
        std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::vector<int> firstRes = firstStep(v, width, height);
  std::vector<int> res = secondStep(firstRes, width, height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      std::cout << res[i * width + j];
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  ASSERT_EQ(right, res);
}
TEST(Sequential_Test, Convex_Test) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  false, true,  true,  true, true,
                      false, true,  false, false, false, false, true, true,
                      true,  true,  false, true,  false, false, true, false,
                      false, true,  true,  true,  false, true};
  std::vector<int> right{135, 0, 135, 0, 137, 135, 2, 135, 0, 3,
                         0,   0, 0,   0, 3,   136, 4, 136, 0, 137,
                         0,   0, 4,   0, 0,   136, 4, 136, 0, 138};
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (v[i * width + j] == true)
        std::cout << static_cast<char>(219);
      else
        std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::vector<int> firstRes = firstStep(v, width, height);
  std::vector<int> secondRes = secondStep(firstRes, width, height);
  std::vector<int> res = convexHull(secondRes, width, height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (res[i * width + j] > 10) {
        std::cout << static_cast<char>(res[i * width + j] - 100);
      } else {
        std::cout << res[i * width + j];
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  ASSERT_EQ(right, res);
}
TEST(Sequential_Test, Convex_Quick_Hull_Test) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  false, true,  true,  true, true,
                      false, true,  false, false, false, false, true, true,
                      true,  true,  false, true,  false, false, true, false,
                      false, true,  true,  true,  false, true};
  std::vector<int> right{139, 0, 139, 0, 140, 139, 2, 139, 0, 3,
                         0,   0, 0,   0, 3,   141, 4, 141, 0, 140,
                         0,   0, 4,   0, 0,   141, 4, 141, 0, 142};
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (v[i * width + j] == true)
        std::cout << static_cast<char>(219);
      else
        std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::vector<int> firstRes = firstStep(v, width, height);
  std::vector<int> secondRes = secondStep(firstRes, width, height);
  std::vector<int> res = convexQuickHull(secondRes, width, height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (res[i * width + j] > 10) {
        std::cout << static_cast<char>(res[i * width + j] - 100);
      } else {
        std::cout << res[i * width + j];
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  ASSERT_EQ(right, res);
}
TEST(Parallel_Test, Convex_Quick_Hull_Test) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  false, true,  true,  true, true,
                      false, true,  false, false, false, false, true, true,
                      true,  true,  false, true,  false, false, true, false,
                      false, true,  true,  true,  false, true};
  std::vector<int> right{143, 0, 143, 0, 144, 143, 2, 143, 0, 3,
                         0,   0, 0,   0, 3,   145, 4, 145, 0, 144,
                         0,   0, 4,   0, 0,   145, 4, 145, 0, 146};
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (v[i * width + j] == true)
        std::cout << static_cast<char>(219);
      else
        std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::vector<int> firstRes = firstStep(v, width, height);
  std::vector<int> secondRes = secondStep(firstRes, width, height);
  std::vector<int> res = convexQuickHullTBB(secondRes, width, height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (res[i * width + j] > 10) {
        std::cout << static_cast<char>(res[i * width + j] - 100);
      } else {
        std::cout << res[i * width + j];
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  ASSERT_EQ(right, res);
}
TEST(Parallel_Test, Convex_Quick_Hull_Test_2) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  true,  true, false, false, false,
                      false, true,  true,  false, true, true,  true,  true,
                      false, false, false, false, true, false, true,  false,
                      true,  true,  false, true,  true, true};
  std::vector<int> right{147, 0, 148, 2, 148, 0,   0, 0,   0, 2,
                         149, 0, 148, 2, 148, 3,   0, 0,   0, 0,
                         3,   0, 150, 0, 150, 149, 0, 150, 6, 150};
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (v[i * width + j] == true)
        std::cout << static_cast<char>(219);
      else
        std::cout << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::vector<int> firstRes = firstStep(v, width, height);
  std::vector<int> secondRes = secondStep(firstRes, width, height);
  std::vector<int> res = convexQuickHullTBB(secondRes, width, height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      if (res[i * width + j] > 10) {
        std::cout << static_cast<char>(res[i * width + j] - 80);
      } else {
        std::cout << res[i * width + j];
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  ASSERT_EQ(right, res);
}

// int main() {
// int width = 500, height = 500;

// std::vector<bool> bitmap = GenerateConnectedMap(width, height);
// /*for (int i = 0; i < height; i++) {
//   for (int j = 0; j < width; j++) {
//     if (bitmap[i * width + j] == true)
//       std::cout << (char)219;
//     else
//       std::cout << " ";
//   }
//   std::cout << "\n";
// }
// std::cout << "\n";*/

//   tbb::tick_count t0 = tbb::tick_count::now();
//   std::vector<int> first1 = firstStep(bitmap, width, height);
//   std::vector<int> second1 = secondStep(first1, width, height);
//   tbb::tick_count t1 = tbb::tick_count::now();
//   std::vector<int> third1 = convexQuickHull(second1, width, height);
//   tbb::tick_count t2 = tbb::tick_count::now();
//   std::vector<int> third2 = convexQuickHullTBB(second1, width, height);
//   tbb::tick_count t3 = tbb::tick_count::now();
//
//   std::cout << "Marking: " << (t1 - t0).seconds() << "\n";
//   std::cout << "QuickConvex: " << (t2 - t1).seconds() << "\n";
//   std::cout << "QuickConvexTBB: " << (t3 - t2).seconds() << "\n\n";
//
//   /*std::cout << "Seq first:\n";
//   for (int i = 0; i < height; i++) {
//     for (int j = 0; j < width; j++) {
//       std::cout << first1[i * width + j];
//     }
//     std::cout << "\n";
//   }
//   std::cout << "\n";
//   std::cout << "Seq second:\n";
//   for (int i = 0; i < height; i++) {
//     for (int j = 0; j < width; j++) {
//       std::cout << second1[i * width + j];
//     }
//     std::cout << "\n";
//   }
//   std::cout << "\n";
//   std::cout << "Seq conv:\n";
//   for (int i = 0; i < height; i++) {
//     for (int j = 0; j < width; j++) {
//       if (third1[i * width + j] > 10)
//         std::cout << (char)(third1[i * width + j] - 100);
//       else
//         std::cout << third1[i * width + j];
//     }
//     std::cout << "\n";
//   }
//   std::cout << "\n";
//   std::cout << "Par conv:\n";
//   for (int i = 0; i < height; i++) {
//     for (int j = 0; j < width; j++) {
//       if (third2[i * width + j] > 10)
//         std::cout << (char)(third2[i * width + j] - 100);
//       else
//         std::cout << third2[i * width + j];
//     }
//     std::cout << "\n";
//   }*/

// return 0;
// }
