// Copyright 2021 Oganyan Robert


#include <gtest/gtest.h>
#include <random>
#include "../../modules/task_4/oganyan_r_mark_components_std/mark_components_std.h"
#include "../../3rdparty/unapproved/unapproved.h"

//  static int count = 0;  //  amount of tests

void Create_Custom_Test(int height, int width) {
    //  ++count;
    //  std::cout<<count<<std::endl;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<unsigned int> distribution(0, 1);
    std::vector<int> seq(width * height);
    for (int i = 0; i < width * height; ++i) {
        seq[i] = distribution(generator);
    }
    auto paral = seq;

    std::cout << height << "x" << width << ": \n";

    auto start = std::chrono::high_resolution_clock::now();
    auto ans = MarkComponentsSeq(&seq, height, width);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "Sequential: " << std::setw(9) << diff.count() << "\n";
    std::cout << " \n";

    start = std::chrono::high_resolution_clock::now();
    auto ans2 = MarkComponentsParStd(&paral, height, width, 4);
    end = std::chrono::high_resolution_clock::now();

    diff = end - start;
    std::cout << "Parallel: " << std::setw(9) << diff.count() << "\n";

    std::cout << " \n";

    ASSERT_EQ(ans.second, ans2.second);
}

TEST(Comparing_to_a_sample, rectangle_image) {
    int width = 3;
    int height = 5;
    std::vector<int> source_image{
            0, 1, 1,
            0, 0, 0,
            1, 1, 1,
            0, 1, 0,
            0, 0, 0,
    };
    std::vector<int> result_image{
            0, 2, 2,
            0, 0, 0,
            3, 3, 3,
            0, 3, 0,
            0, 0, 0,
    };
    std::pair<std::vector<int>, int> res_image_seq;
    std::pair<std::vector<int>, int> res_image_paral;
    ASSERT_NO_THROW(res_image_seq = MarkComponentsSeq(&source_image, height, width));
    ASSERT_NO_THROW(res_image_paral = MarkComponentsParStd(&source_image, height, width, 4));
    ASSERT_EQ(res_image_seq.second, 2);
    ASSERT_EQ(res_image_seq.second, res_image_paral.second);
}

TEST(Wrong_data, wrong_img_2) {
    std::vector<int> source_image{};
    std::vector<int> result_image{};
    int width = 5;
    int height = 2;
    ASSERT_ANY_THROW(MarkComponentsSeq(&source_image, height, width));
    ASSERT_ANY_THROW(MarkComponentsParStd(&source_image, height, width, 4));
}


TEST(Comparing_to_a_sample, square_image_1) {
    int width = 6;
    int height = 6;
    std::vector<int> source_image{
            0, 1, 0, 0, 0, 0,
            1, 1, 1, 0, 0, 1,
            0, 0, 0, 0, 1, 1,
            0, 1, 0, 1, 0, 0,
            1, 1, 1, 0, 0, 0,
            0, 0, 0, 0, 1, 1,
    };
    std::vector<int> result_image{
            0, 2, 0, 0, 0, 0,
            2, 2, 2, 0, 0, 3,
            0, 0, 0, 0, 3, 3,
            0, 4, 0, 5, 0, 0,
            4, 4, 4, 0, 0, 0,
            0, 0, 0, 0, 6, 6,
    };
    std::pair<std::vector<int>, int> res_image_seq;
    std::pair<std::vector<int>, int> res_image_paral;
    ASSERT_NO_THROW(res_image_seq = MarkComponentsSeq(&source_image, height, width));
    ASSERT_NO_THROW(res_image_paral = MarkComponentsParStd(&source_image, height, width, 4));
    ASSERT_EQ(res_image_seq.second, 5);
    ASSERT_EQ(res_image_seq.second, res_image_paral.second);
}

TEST(Wrong_data, wrong_width_and_height_3) {
    std::vector<int> source_image{
            0, 0, 1,
            1, 0, 1,
            0, 0, 0,
    };
    std::vector<int> result_image{
            0, 0, 2,
            3, 0, 2,
            0, 0, 0,
    };
    int width = 1;
    int height = 3;
    ASSERT_ANY_THROW(MarkComponentsSeq(&source_image, height, width));
    ASSERT_ANY_THROW(MarkComponentsParStd(&source_image, height, width, 4));
}

TEST(Comparing_to_a_big_sample, square_image_1) {
    int width = 15;
    int height = 15;
    std::vector<int> source_image{
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    };
    std::vector<int> result_image{
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    };
    std::pair<std::vector<int>, int> res_image_seq;
    std::pair<std::vector<int>, int> res_image_paral;
    ASSERT_NO_THROW(res_image_seq = MarkComponentsSeq(&source_image, height, width));
    ASSERT_NO_THROW(res_image_paral = MarkComponentsParStd(&source_image, height, width, 4));
    ASSERT_EQ(res_image_seq.second, 1);
    ASSERT_EQ(res_image_seq.second, res_image_paral.second);
}

TEST(Wrong_data, wrong_width_and_height_1) {
    std::vector<int> source_image = {
            0, 1, 1,
            1, 0, 0,
            0, 0, 0,
    };
    std::vector<int> result_image = {
            0, 2, 2,
            3, 0, 0,
            0, 0, 0,
    };
    int width = 4;
    int height = 3;
    ASSERT_ANY_THROW(MarkComponentsSeq(&source_image, height, width));
    ASSERT_ANY_THROW(MarkComponentsParStd(&source_image, height, width, 4));
}

TEST(Comparing_to_a_big_sample, square_image_2) {
    int width = 10;
    int height = 10;
    std::vector<int> source_image{
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
            1, 0, 1, 0, 0, 0, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 0, 1,
            1, 0, 1, 0, 1, 1, 0, 1, 0, 1,
            1, 0, 1, 0, 0, 0, 0, 1, 0, 1,
            1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    };
    std::vector<int> result_image{
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 0, 0, 0, 0, 0, 0, 0, 0, 2,
            2, 0, 3, 3, 3, 3, 3, 3, 0, 2,
            2, 0, 3, 0, 0, 0, 0, 3, 0, 2,
            2, 0, 3, 0, 4, 4, 0, 3, 0, 2,
            2, 0, 3, 0, 4, 4, 0, 3, 0, 2,
            2, 0, 3, 0, 0, 0, 0, 3, 0, 2,
            2, 0, 3, 3, 3, 3, 3, 3, 0, 2,
            2, 0, 0, 0, 0, 0, 0, 0, 0, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    };
    std::pair<std::vector<int>, int> res_image_seq;
    std::pair<std::vector<int>, int> res_image_paral;
    ASSERT_NO_THROW(res_image_seq = MarkComponentsSeq(&source_image, height, width));
    ASSERT_NO_THROW(res_image_paral = MarkComponentsParStd(&source_image, height, width, 4));
    ASSERT_EQ(res_image_seq.second, 3);
    ASSERT_EQ(res_image_seq.second, res_image_paral.second);
}

TEST(Wrong_data, wrong_img) {
    std::vector<int> source_image;
    std::vector<int> result_image;
    int width = 0;
    int height = 0;
    ASSERT_ANY_THROW(MarkComponentsSeq(&source_image, height, width));
    ASSERT_ANY_THROW(MarkComponentsParStd(&source_image, height, width, 4));
}


TEST(Comparing_to_a_sample, square_image_2) {
    int width = 5;
    int height = 5;
    std::vector<int> source_image{
            0, 1, 0, 1, 1,
            0, 0, 0, 1, 1,
            1, 0, 0, 0, 1,
            0, 0, 1, 1, 1,
            0, 0, 0, 0, 1,
    };
    std::vector<int> result_image{
            0, 2, 0, 3, 3,
            0, 0, 0, 3, 3,
            4, 0, 0, 0, 3,
            0, 0, 3, 3, 3,
            0, 0, 0, 0, 3,
    };
    std::pair<std::vector<int>, int> res_image_seq;
    std::pair<std::vector<int>, int> res_image_paral;
    ASSERT_NO_THROW(res_image_seq = MarkComponentsSeq(&source_image, height, width));
    ASSERT_NO_THROW(res_image_paral = MarkComponentsParStd(&source_image, height, width, 4));
    ASSERT_EQ(res_image_seq.second, 3);
    ASSERT_EQ(res_image_seq.second, res_image_paral.second);
}

TEST(Wrong_data, wrong_img_1) {
    std::vector<int> source_image{};
    std::vector<int> result_image{};
    int width = 0;
    int height = 0;
    ASSERT_ANY_THROW(MarkComponentsSeq(&source_image, height, width));
    ASSERT_ANY_THROW(MarkComponentsParStd(&source_image, height, width, 4));
}

TEST(Wrong_data, wrong_width_and_height_2) {
    std::vector<int> source_image{
            0, 1, 1,
            0, 0, 1,
            1, 1,
    };
    std::vector<int> result_image{
            0, 2, 2,
            0, 0, 2,
            2, 2,
    };
    int width = 3;
    int height = 3;
    ASSERT_ANY_THROW(MarkComponentsSeq(&source_image, height, width));
    ASSERT_ANY_THROW(MarkComponentsParStd(&source_image, height, width, 4));
}


TEST(Comparing_to_a_big_sample, square_image_3) {
    int width = 10;
    int height = 10;
    std::vector<int> source_image{
            0, 1, 0, 1, 0, 1, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
            0, 0, 1, 1, 0, 1, 0, 0, 1, 0,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 1, 0, 1, 1, 1, 1, 0, 0,
            1, 0, 1, 1, 1, 1, 1, 0, 0, 1,
            0, 1, 1, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 1, 1, 0, 1, 0, 1, 1,
            1, 0, 1, 0, 0, 0, 0, 0, 1, 1,
            1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
    };
    std::vector<int> result_image{
            0, 2, 0, 3, 0, 4, 4, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 4, 4, 0, 0,
            0, 0, 5, 5, 0, 6, 0, 0, 7, 0,
            8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            8, 0, 9, 0, 9, 9, 9, 9, 0, 0,
            8, 0, 9, 9, 9, 9, 9, 0, 0, 10,
            0, 9, 9, 0, 0, 0, 0, 11, 0, 0,
            0, 0, 0, 12, 12, 0, 13, 0, 14, 14,
            15, 0, 16, 0, 0, 0, 0, 0, 14, 14,
            15, 0, 16, 0, 17, 0, 18, 0, 14, 14,
    };

    std::pair<std::vector<int>, int> res_image_seq;
    std::pair<std::vector<int>, int> res_image_paral;
    ASSERT_NO_THROW(res_image_seq = MarkComponentsSeq(&source_image, height, width));
    ASSERT_NO_THROW(res_image_paral = MarkComponentsParStd(&source_image, height, width, 4));
    ASSERT_EQ(res_image_seq.second, 17);
    ASSERT_EQ(res_image_seq.second, res_image_paral.second);
}


TEST(Comparing_random_img, small_image) {
    Create_Custom_Test(10, 10);
}

TEST(Comparing_random_img, small_image_2) {
    Create_Custom_Test(100, 100);
}

TEST(Comparing_random_img, middle_image) {
    Create_Custom_Test(500, 500);
}

TEST(Comparing_random_img, middle_image_2) {
    Create_Custom_Test(750, 750);
}

TEST(Comparing_random_img, middle_image_3) {
    Create_Custom_Test(1000, 1000);
}

/*

TEST(Comparing_random_img, big_image) {
    Create_Custom_Test(2000, 2000);
}

TEST(Comparing_random_img, big_image_2) {
    Create_Custom_Test(3000, 3000);
}

TEST(Comparing_random_img, big_image_3) {
    Create_Custom_Test(5000, 5000);
}

*/

int main() {
    testing::InitGoogleTest();
//    testing::FLAGS_gtest_repeat = 40;
    return RUN_ALL_TESTS();
}
