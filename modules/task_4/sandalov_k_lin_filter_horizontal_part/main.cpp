// Copyright 2020 Sandalov Konstantin
#include <gtest/gtest.h>
#include <vector>
#include "../../modules/task_4/sandalov_k_lin_filter_horizontal_part/lin_filter_horizontal_part.h"


static uint8_t init_arr[] {
    255, 5, 255, 254, 255,
    255, 255, 252, 255, 255,
    255, 254, 1, 252, 253,
    253, 255, 252, 255, 0,
    254, 3, 253, 255, 254
};
static my::Image<uint8_t> testImage(init_arr, 5, 5, 1);

static uint8_t blank_image[] {
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
};
static my::Image<uint8_t> blankImage(blank_image, 5, 5, 1);

static const std::vector<uint8_t> res_arr {
    205, 172, 204, 254, 254,
    236, 204, 203, 234, 254,
    254, 222, 201, 203, 203,
    235, 203, 203, 203, 170,
    203, 171, 203, 235, 203
};

static const std::vector<float> correct_kernel {
    0.0751136, 0.123841, 0.0751136,
    0.123841, 0.20418, 0.123841,
    0.0751136, 0.123841, 0.0751136
};

TEST(STD_Filter, Create_Filter_With_Even_Size) {
    ASSERT_ANY_THROW(my::createGaussianFilter(6, 1));
}

TEST(STD_Filter, Create_Filter_With_Negative_Sigma) {
    ASSERT_ANY_THROW(my::createGaussianFilter(3, -1));
}

TEST(STD_Filter, Create_Filter_With_Correct_Params) {
    std::vector<float> computedKernel;
    ASSERT_NO_THROW(computedKernel = my::createGaussianFilter(3, 1));
    ASSERT_EQ(computedKernel.size(), correct_kernel.size());
    for (size_t i = 0; i < computedKernel.size(); ++i) {
        ASSERT_NEAR(computedKernel[i], correct_kernel[i], 0.000001);
    }
}

TEST(STD_Filter, Image_Processing_With_Incorrect_Kernel_Size) {
    std::vector<float> filter(9, 0);
    ASSERT_ANY_THROW(my::linearFilterSeq(testImage, filter, 5));
    ASSERT_ANY_THROW(my::linearFilterParStd(testImage, filter, 5));
}

TEST(STD_Filter, Filter_Doesnt_Change_Blank_Image) {
    auto kernel = my::createGaussianFilter(3, 1);
    auto resSeq = my::linearFilterSeq(blankImage, kernel, 3);
    for (size_t i = 0; i < static_cast<size_t>(blankImage.width*blankImage.height*blankImage.dims); ++i) {
        ASSERT_EQ(blank_image[i], resSeq[i]);
    }
    auto resPar = my::linearFilterParStd(blankImage, kernel, 3);
    for (size_t i = 0; i < static_cast<size_t>(blankImage.width*blankImage.height*blankImage.dims); ++i) {
        ASSERT_EQ(blank_image[i], resPar[i]);
    }
    ASSERT_EQ(resSeq.size(), resPar.size());
    for (size_t i = 0; i < resSeq.size(); ++i) {
        ASSERT_EQ(resSeq[i], resPar[i]);
    }
}

TEST(STD_Filter, Filetr_Work_Correct) {
    auto kernel = my::createGaussianFilter(3, 1);
    auto resSeq = my::linearFilterSeq(testImage, kernel, 3);
    for (size_t i = 0; i < resSeq.size(); ++i) {
        ASSERT_EQ(res_arr[i], resSeq[i]);
    }
    auto resPar = my::linearFilterParStd(testImage, kernel, 3);
    for (size_t i = 0; i < resPar.size(); ++i) {
        ASSERT_EQ(res_arr[i], resPar[i]);
    }
    ASSERT_EQ(resSeq.size(), resPar.size());
    for (size_t i = 0; i < resSeq.size(); ++i) {
        ASSERT_EQ(resSeq[i], resPar[i]);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
