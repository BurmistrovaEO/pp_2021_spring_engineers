// Copyright 2021 Ognev Denis

#include <gtest/gtest.h>

#include "../../task_2/ognev_d_markings_hulls_omp/markings_hulls_omp.h"

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
TEST(Parallel_Test, Marking_Test) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  false, true,  true,  true, true,
                      false, true,  false, false, false, false, true, true,
                      true,  true,  false, true,  false, false, true, false,
                      false, true,  true,  true,  false, true};
  std::vector<int> right{2, 0, 2, 0, 4, 2, 2, 2, 0, 4, 0, 0, 0, 0, 4,
                         5, 5, 5, 0, 4, 0, 0, 5, 0, 0, 5, 5, 5, 0, 6};
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
  std::vector<int> firstRes = firstStepOMPver2(v, width, height);
  std::vector<int> secondRes = secondStepOMPver2(firstRes, width, height);
  ASSERT_NO_THROW(secondStepOMPver2(firstStepOMPver2(v, width, height), width, height));
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      std::cout << secondRes[i * width + j];
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}
TEST(Parallel_Test, Convex_QuickHull_Test) {
  int width = 5, height = 6;
  std::vector<bool> v{true,  false, true,  false, true,  true,  true, true,
                      false, true,  false, false, false, false, true, true,
                      true,  true,  false, true,  false, false, true, false,
                      false, true,  true,  true,  false, true};
  std::vector<int> right{143, 0, 143, 0, 144, 143, 2, 143, 0, 4,
                         0,   0, 0,   0, 4,   145, 5, 145, 0, 144,
                         0,   0, 5,   0, 0,   145, 5, 145, 0, 146};
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
  std::vector<int> firstRes = firstStepOMPver2(v, width, height);
  std::vector<int> secondRes = secondStepOMPver2(firstRes, width, height);
  std::vector<int> res = convexQuickHullOMP(secondRes, width, height);
  ASSERT_NO_THROW(convexQuickHullOMP(
      secondStepOMPver2(firstStepOMP(v, width, height), width, height), width,
      height));

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
}

// int main() {
//  std::vector<int> treeUsual(1, 0);
//  std::vector<int> treeUsual2(1, 0);
//  std::vector<int> treeOMP(1, 0);
//
//  int width = 5, height = 6;
//
//  std::vector<bool> bitmap = GenerateConnectedMap(width, height);
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      if (bitmap[i * width + j] == true)
//        std::cout << (char)219;
//      else
//        std::cout << " ";
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//
//  // omp_set_nested(1);
//
//  double t1 = omp_get_wtime();
//  std::vector<int> first1 = firstStep(bitmap, width, height, treeUsual);
//  std::vector<int> second1 = secondStep(first1, width, height, treeUsual);
//  double t2 = omp_get_wtime();
//  std::vector<int> third1 = convexQuickHull(second1, width, height);
//  double t3 = omp_get_wtime();
//  std::vector<int> first2 = firstStep(bitmap, width, height, treeOMP);
//  std::vector<int> second2 = secondStep(first2, width, height, treeOMP);
//  double t4 = omp_get_wtime();
//  std::vector<int> third2 = convexQuickHullOMP(second2, width, height);
//  double t5 = omp_get_wtime();
//  std::vector<int> first3 = firstStepOMP(bitmap, width, height, treeUsual2);
//  std::vector<int> second3 = secondStepOMP(first3, width, height, treeUsual2);
//  double t6 = omp_get_wtime();
//
//  std::cout << "Marking: " << t2 - t1 << "\n";
//  std::cout << "MarkingOMP: " << t6 - t5 << "\n";
//  std::cout << "QuickConvex: " << t3 - t2 << "\n";
//  std::cout << "QuickConvexOMP: " << t5 - t4 << "\n\n";
//  std::cout << "Sequential time: " << t3 - t1 << "\n";
//  std::cout << "Parallel time: " << t5 - t3 << "\n\n";
//
//  std::cout << "Seq first:\n";
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      std::cout << first1[i * width + j];
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//  std::cout << "Par first:\n";
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      std::cout << first3[i * width + j];
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//  std::cout << "Seq second:\n";
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      std::cout << second1[i * width + j];
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//  std::cout << "Par second:\n";
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      std::cout << second2[i * width + j];
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//  std::cout << "Seq conv:\n";
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      if (third1[i * width + j] > 10)
//        std::cout << (char)(third1[i * width + j] - 100);
//      else
//        std::cout << third1[i * width + j];
//    }
//    std::cout << "\n";
//  }
//  std::cout << "\n";
//  std::cout << "Par conv:\n";
//  for (int i = 0; i < height; i++) {
//    for (int j = 0; j < width; j++) {
//      if (third2[i * width + j] > 10)
//        std::cout << (char)(third2[i * width + j] - 100);
//      else
//        std::cout << third2[i * width + j];
//    }
//    std::cout << "\n";
//  }
//  return 0;
// }
