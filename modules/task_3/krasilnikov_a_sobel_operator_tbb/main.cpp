// Copyright 2021 Krasilnikov Alexey
#include <gtest/gtest.h>
#include <random>
#include "../../modules/task_3/krasilnikov_a_sobel_operator_tbb/sobel_operator_tbb.h"


TEST(SobelOperator, Smooth) {
    std::vector<unsigned char> img = {0, 0, 0,
                                      0, 0, 0,
                                      0, 0, 0};
    SobelOperator(img.data(), 3, 3);
    ASSERT_EQ(9u, img.size());
    ASSERT_EQ(0, img[4]);

    img = {255, 255, 255,
           255, 255, 255,
           255, 255, 255};
    SobelOperator(img.data(), 3, 3);
    ASSERT_EQ(9u, img.size());
    ASSERT_EQ(0, img[4]);
}

TEST(SobelOperator, MaximumDifference) {
    std::vector<unsigned char> img = {255, 0, 0,
                                      255, 0, 0,
                                      255, 0, 0};
    SobelOperator(img.data(), 3, 3);
    ASSERT_EQ(9u, img.size());
    ASSERT_EQ(255, img[4]);

    img = {255, 255, 255,
           0, 0, 0,
           0, 0, 0};
    SobelOperator(img.data(), 3, 3);
    ASSERT_EQ(9u, img.size());
    ASSERT_EQ(255, img[4]);

    img = {0, 0, 0,
           0, 255, 0,
           255, 255, 255};
    SobelOperator(img.data(), 3, 3);
    ASSERT_EQ(9u, img.size());
    ASSERT_EQ(255, img[4]);
}

TEST(SobelOperator, Arithmetic) {
    std::vector<unsigned char> img = {0, 10, 20,
                                      70, 0, 30,
                                      60, 50, 40};
    SobelOperator(img.data(), 3, 3);
    ASSERT_EQ(9u, img.size());
    ASSERT_EQ(179, img[4]);
}

TEST(SobelOperator, parallel_1) {
    std::random_device rd;
    std::mt19937 mersenne(rd());

    std::vector<unsigned char> img(100*100);
    std::vector<unsigned char> img2(100*100);
    for (int i = 0; i < 100 * 100; ++i) {
        img[i] = mersenne() % 256;
        img2[i] = img[i];
    }
    SobelOperator(img.data(), 100, 100);
    SobelOperator_parallel(img2.data(), 100, 100);
    for (int i = 0; i < 100 * 100; ++i) {
        ASSERT_EQ(img[i], img2[i]);
    }
}

TEST(SobelOperator, parallel_2) {
    std::random_device rd;
    std::mt19937 mersenne(rd());

    std::vector<unsigned char> img(1000*1000);
    std::vector<unsigned char> img2(1000*1000);
    for (int i = 0; i < 1000 * 1000; ++i) {
        img[i] = mersenne() % 256;
        img2[i] = img[i];
    }
    SobelOperator(img.data(), 1000, 1000);
    SobelOperator_parallel(img2.data(), 1000, 1000);
    for (int i = 0; i < 1000 * 1000; ++i) {
        ASSERT_EQ(img[i], img2[i]);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
