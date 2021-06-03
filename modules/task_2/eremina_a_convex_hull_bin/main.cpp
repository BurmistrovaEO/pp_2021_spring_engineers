// Copyright 2021 Eremina Alena
#include <gtest/gtest.h>
#include <ctime>
#include "./convex_hull_bin.h"

TEST(Convex_Hull_Bin, Test_OMP_Convex_Bin_Random_Image_No_Throw) {
    int height = 20;
    int width = 40;
    std::vector<uint8_t> img_src(height * width);
    std::vector<uint8_t> img_dst(height * width, 255);
    img_src = generateImage(width, height);
    ASSERT_NO_THROW(img_dst = Convex_Hull_OMP(img_src, height, width));
}

TEST(Convex_Hull_Bin, Test_OMP_Convex_Bin_Random_Image) {
    int height = 20;
    int width = 40;
    std::vector<uint8_t> img_src(height * width);
    std::vector<uint8_t> img_dst_seq(height * width, 255);
    std::vector<uint8_t> img_dst_omp(height * width, 255);
    img_src = generateImage(width, height);
    img_dst_seq = Convex_Hull(img_src, height, width);
    img_dst_omp = Convex_Hull_OMP(img_src, height, width);
    int count = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (img_dst_seq[i * width + j] == img_dst_omp[i * width + j]) {
                count++;
            }
        }
    }
    double ratio = static_cast<double>(count)/(width*height);
    ASSERT_TRUE(ratio >= 0.95);
}

TEST(Convex_Hull_Bin, Test_Convex_Bin_Small_Image) {
    int height = 8;
    int width = 8;
    std::vector<uint8_t> img_src = { 255, 255, 255, 255, 255, 255, 255, 255,
                                     255, 255, 255,   0, 255, 255, 255, 255,
                                     255, 255, 255,   0, 255, 255, 255, 255,
                                     255, 255, 255,   0,   0, 255, 255, 255,
                                     255, 255,   0,   0,   0,   0, 255, 255,
                                     255,   0,   0,   0,   0,   0,   0, 255,
                                     255, 255,   0,   0,   0,   0, 255, 255,
                                     255, 255, 255, 255, 255, 255, 255, 255 };
    std::vector<uint8_t> img_dst_seq(height * width);
    std::vector<uint8_t> img_dst_omp(height * width);
    img_dst_seq = Convex_Hull(img_src, height, width);
    img_dst_omp = Convex_Hull_OMP(img_src, height, width);
    int count = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (img_dst_seq[i * width + j] == img_dst_omp[i * width + j]) {
                count++;
            }
        }
    }
    double ratio = static_cast<double>(count)/(width*height);
    ASSERT_TRUE(ratio >= 0.95);
}

TEST(Convex_Hull_Bin, Test_Convex_Bin_Medium_Image) {
    int height = 10;
    int width = 20;
    std::vector<uint8_t> img_src =
        { 255, 255, 255, 255,   0,   0, 255, 255, 255, 255,
          255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
          255, 255, 255, 255,   0,   0, 255, 255, 255, 255,
          255,   0, 255, 255, 255, 255, 255, 255,   0, 255,
          255, 255, 255,   0,   0,   0, 255, 255, 255, 255,
          255, 255,   0, 255, 255, 255, 255,   0, 255, 255,
          255, 255, 255,   0,   0,   0, 255, 255, 255, 255,
          255, 255, 255,   0, 255, 255,   0, 255, 255, 255,
            0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
          255, 255, 255, 255,   0,   0, 255, 255, 255, 255,
            0,   0,   0,   0,   0,   0,   0,   0, 255, 255,
          255, 255, 255, 255,   0,   0, 255, 255, 255, 255,
          255, 255, 255,   0,   0,   0,   0, 255, 255, 255,
          255, 255, 255,   0, 255, 255,   0, 255, 255, 255,
          255, 255, 255, 255,   0,   0,   0, 255, 255, 255,
          255, 255,   0, 255, 255, 255, 255,   0, 255, 255,
          255, 255, 255, 255,   0, 255, 255, 255, 255, 255,
          255,   0, 255, 255, 255, 255, 255, 255,   0, 255,
          255, 255, 255, 255,   0, 255, 255, 255, 255, 255,
          255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };
    std::vector<uint8_t> img_dst_seq(height * width);
    std::vector<uint8_t> img_dst_omp(height * width);
    img_dst_seq = Convex_Hull(img_src, height, width);
    img_dst_omp = Convex_Hull_OMP(img_src, height, width);
    int count = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (img_dst_seq[i * width + j] == img_dst_omp[i * width + j]) {
                count++;
            }
        }
    }
    double ratio = static_cast<double>(count)/(width*height);
    ASSERT_TRUE(ratio >= 0.95);
}

TEST(Convex_Hull_Bin, Test_Convex_Bin_Large_Image) {
    int height = 32;
    int width = 20;
    std::vector<uint8_t> img_src =
    { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255,   0,   0,   0,   0,   0,   0,   0,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255,   0,   0,   0,
      255, 255,   0,   0,   0,   0, 255, 255, 255, 255,
        0, 255, 255, 255, 255, 255, 255, 255,   0,   0,
      255, 255,   0,   0,   0,   0,   0,   0,   0, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255,   0,
      255,   0,   0,   0,   0,   0,   0,   0,   0, 255,
      255,   0, 255, 255, 255, 255, 255, 255, 255,   0,
      255,   0,   0,   0,   0,   0,   0,   0,   0, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255,   0,
      255,   0,   0,   0, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0,   0, 255, 255, 255, 255, 255, 255,
        0, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0,   0, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0, 255, 255, 255, 255, 255,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      255, 255, 255, 255, 255, 255, 255, 255,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        0,   0,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255,   0,
        0,   0,   0,   0,   0,   0,   0,   0, 255, 255,
      255, 255, 255, 255, 255,   0,   0, 255,   0,   0,
        0,   0,   0,   0,   0,   0,   0, 255, 255, 255,
      255, 255, 255, 255,   0,   0,   0, 255,   0, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255,   0,   0, 255,   0,   0,   0, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255,   0, 255, 255,   0,   0, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255,   0,   0, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255,   0,   0, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

    std::vector<uint8_t> img_dst_seq(height * width);
    std::vector<uint8_t> img_dst_omp(height * width);
    img_dst_seq = Convex_Hull(img_src, height, width);
    img_dst_omp = Convex_Hull_OMP(img_src, height, width);
    int count = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (img_dst_seq[i * width + j] == img_dst_omp[i * width + j]) {
                count++;
            }
        }
    }
    double ratio = static_cast<double>(count)/(width*height);
    ASSERT_TRUE(ratio >= 0.95);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
