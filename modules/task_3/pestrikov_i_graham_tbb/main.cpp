// Copyright 2021 Pestrikov Ilya

#include <gtest/gtest.h>
#include <tbb/tbb.h>
#include "../../modules/task_3/pestrikov_i_graham_tbb/graham.h"

TEST(Graham, Test_Square) {
    std::vector<Point> points = std::vector<Point>();
    points.push_back({ 2, 2 });
    points.push_back({ 2, 1 });
    points.push_back({ 1, 2 });
    points.push_back({ -1, -1 });

    std::vector<Point>* grahamResult = graham(&points, 0, points.size());
    std::vector<Point> expectedResult = std::vector<Point>();
    expectedResult.push_back({ -1, -1 });
    expectedResult.push_back({ 2, 1 });
    expectedResult.push_back({ 2, 2 });
    expectedResult.push_back({ 1, 2 });

    ASSERT_EQ(grahamResult->size(), expectedResult.size());

    for (size_t i = 0; i < grahamResult->size(); ++i) {
        ASSERT_EQ(grahamResult->operator[](i).x, expectedResult[i].x);
        ASSERT_EQ(grahamResult->operator[](i).y, expectedResult[i].y);
    }
}

TEST(Graham, Test_Random_10_Points_In_Range_10_20) {
    auto points = generateRandomPointsVector(10, 10, 20);

    clock_t start = clock();
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    clock_t end = clock();
    std::cout << "Graham sequential time: " <<
        static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000 <<
        " ms" << std::endl;

    std::vector<Point>* grahamTbbResult =
        grahamTbb(new std::vector<Point>(*points));
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamTbbResult, *points));
}

TEST(Graham, Test_Random_100_Points_In_Range_100_200) {
    auto points = generateRandomPointsVector(100, 100, 200);

    clock_t start = clock();
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    clock_t end = clock();
    std::cout << "Graham sequential time: " <<
        static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000 <<
        " ms" << std::endl;

    std::vector<Point>* grahamTbbResult =
        grahamTbb(new std::vector<Point>(*points));
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamTbbResult, *points));
}

TEST(Graham, Test_Random_1000_Points_In_Range_1000_2000) {
    auto points = generateRandomPointsVector(1000, 1000, 2000);

    clock_t start = clock();
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    clock_t end = clock();
    std::cout << "Graham sequential time: " <<
        static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000 <<
        " ms" << std::endl;

    std::vector<Point>* grahamTbbResult =
        grahamTbb(new std::vector<Point>(*points));
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamTbbResult, *points));
}

TEST(Graham, Test_Random_1_000_000_Points_In_Range_10000_20000) {
    auto points = generateRandomPointsVector(1000000, 10000, 20000);

    clock_t start = clock();
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    clock_t end = clock();
    std::cout << "Graham sequential time: " <<
        static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000 <<
        " ms" << std::endl;

    std::vector<Point>* grahamTbbResult =
        grahamTbb(new std::vector<Point>(*points));
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamTbbResult, *points));
}

// This test takes a lot of time to execute (~40 seconds),
// so should be used only for demonstration purposes

// TEST(Graham, Test_Random_10_000_000_Points_In_Range_100000_200000) {
//    auto points = generateRandomPointsVector(10000000, 100000, 200000);
//
//    clock_t start = clock();
//    std::vector<Point>* grahamResult =
//        graham(new std::vector<Point>(*points), 0, points->size());
//    clock_t end = clock();
//    std::cout << "Graham sequential time: " <<
//        static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000 <<
//        " ms" << std::endl;
//
//    std::vector<Point>* grahamTbbResult =
//        grahamTbb(new std::vector<Point>(*points));
//    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
//    ASSERT_TRUE(isGrahamCorrect(*grahamTbbResult, *points));
//}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
