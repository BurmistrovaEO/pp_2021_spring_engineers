// Copyright 2020 Egorov Kirill
#include <gtest/gtest.h>
#include <omp.h>
#include <vector>
#include "./egorov_k_gauss_filter_omp.h"

TEST(basic_tests, clamp_test) {
    int num[3] = { -1, 2, 3 };

    ASSERT_EQ(clamp(num[0], 0, 2), 0);
    ASSERT_EQ(clamp(num[1], 0, 2), 2);
    ASSERT_EQ(clamp(num[2], 0, 2), 2);
}

TEST(basic_tests, kernel_test) {
    std::vector<float> test(9);

    ASSERT_NO_THROW(test = calculateKernel3x3(3));
}

TEST(basic_tests, parallel_test) {
    std::vector<int> img = getRandomVector(6400);
    std::vector<int> res(img.size());

    ASSERT_NO_THROW(res = filterApplyParallel(img, 80, 80, 3));
}

TEST(basic_tests, seq_test) {
    std::vector<int> img = getRandomVector(6400);
    std::vector<int> res(img.size());

    ASSERT_NO_THROW(res = filterApply(img, 80, 80, 3));
}

TEST(efficiency_tests, omp_seq_test) {
    std::vector<int> img = getRandomVector(6400);
    std::vector<int> res(img.size());
    double s, e, total;
    double sP, eP, totalP;

    s = omp_get_wtime();;
    res = filterApply(img, 80, 80, 3);
    e = omp_get_wtime();;
    total = e - s;

    sP = omp_get_wtime();
    res = filterApplyParallel(img, 80, 80, 3);
    eP = omp_get_wtime();
    totalP = eP - sP;

    std::cout << "Times: p: " << totalP << " s: " << total << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
