// Copyright 2021 Novikov Danil

#include <gtest/gtest.h>
#include <ctime>
#include "../../../modules/task_2/novikov_d_rectangles_method_omp/novikov_d_rectangles_method_omp.h"

TEST(Novikov, TestException1) {
    std::vector <double> begin_point = { 0 };
    std::vector <double> end_point = { 2 };
    int number_of_partitions = 40;

    ASSERT_NO_THROW(rectangles_omp([](std::vector <double> vec){return 0;},
    begin_point, end_point, number_of_partitions));
}

TEST(Novikov, TestException2) {
    std::vector <double> begin_point = { 0 };
    std::vector <double> end_point = { 2 };
    int number_of_partitions = -1;

    ASSERT_ANY_THROW(rectangles_omp([](std::vector <double> vec){return 0;},
    begin_point, end_point, number_of_partitions));
}

TEST(Novikov, TestException3) {
    std::vector <double> begin_point = { 0, 2 };
    std::vector <double> end_point = { 2 };
    int number_of_partitions = 5;

    ASSERT_ANY_THROW(rectangles_omp([](std::vector <double> vec){return 0;},
    begin_point, end_point, number_of_partitions));
}

TEST(Novikov, TestException4) {
    std::vector <double> begin_point = { 0, 0, 0, 0, 0, 0 };
    std::vector <double> end_point = { 2, 2, 2, 2, 2, 2 };
    int number_of_partitions = 5;

    ASSERT_NO_THROW(rectangles_omp([](std::vector <double> vec){return 0;},
    begin_point, end_point, number_of_partitions));
}

TEST(Novikov, TestEazy4d) {
    std::vector <double> begin_point = { 0, 0, 0, 0 };
    std::vector <double> end_point = { 2, 2, 2, 2 };
    double result;
    int number_of_partitions = 20;

    auto function = [](std::vector <double> vec) {
        return vec[0]*vec[0]*3;
    };

    result = rectangles_omp(function, begin_point,
                            end_point, number_of_partitions);
    EXPECT_NEAR(result, 64.0, 0.1);
}

TEST(Novikov, TestEazy3d) {
    std::vector <double> begin_point = { 0, 0, 0 };
    std::vector <double> end_point = { 3, 3, 3 };
    double result;
    int number_of_partitions = 20;

    auto function = [](std::vector <double> vec) {
        return vec[0] + vec[1] + vec[2];
    };

    result = rectangles_omp(function, begin_point,
                            end_point, number_of_partitions);
    EXPECT_NEAR(result, 121.5, 0.5);
}

TEST(Novikov, TestSin2d) {
    std::vector <double> begin_point = { 0, 0 };
    std::vector <double> end_point = { 2, 2 };
    double result;
    int number_of_partitions = 30;

    auto function = [](std::vector <double> vec) {
        return sin(vec[0]) + sin(vec[1]);
    };

    result = rectangles_omp(function, begin_point,
                            end_point, number_of_partitions);
    EXPECT_NEAR(result, 4*(1 - cos(2)), 0.1);
}

TEST(Novikov, TestTimings) {
    std::vector <double> begin_point = { 0, 0, 0 };
    std::vector <double> end_point = { 3, 3, 3 };
    double result1, result2;
    int number_of_partitions = 20;

    auto function = [](std::vector <double> vec) {
        return vec[0] + vec[1] + vec[2];
    };

    clock_t begin_time = clock();
    result1 = rectangles_base(function, begin_point,
                            end_point, number_of_partitions);
    std::cout << "Sequential time: " << static_cast<float>(clock() -
    begin_time) / CLOCKS_PER_SEC << std::endl;

    begin_time = clock();
    result2 = rectangles_omp(function, begin_point,
                            end_point, number_of_partitions);
    std::cout << "OMP time: " << static_cast<float>(clock() -
    begin_time) / CLOCKS_PER_SEC << std::endl;

    EXPECT_NEAR(result1, result2, 0.01);
}
