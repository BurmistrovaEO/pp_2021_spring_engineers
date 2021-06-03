// Copyright 2021 Pestrikov Ilya
#include <gtest/gtest.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include "./pestrikov_i_graham_omp.h"

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
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    std::vector<Point>* grahamOmpResult2Threads =
        grahamOmp(new std::vector<Point>(*points), 2);
    std::vector<Point>* grahamOmpResult3Threads =
        grahamOmp(new std::vector<Point>(*points), 3);
    std::vector<Point>* grahamOmpResult4Threads =
        grahamOmp(new std::vector<Point>(*points), 4);
    std::vector<Point>* grahamOmpResult8Threads =
        grahamOmp(new std::vector<Point>(*points), 8);
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult2Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult3Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult4Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult8Threads, *points));
}

TEST(Graham, Test_Random_100_Points_In_Range_100_200) {
    auto points = generateRandomPointsVector(100, 100, 200);
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    std::vector<Point>* grahamOmpResult2Threads =
        grahamOmp(new std::vector<Point>(*points), 2);
    std::vector<Point>* grahamOmpResult3Threads =
        grahamOmp(new std::vector<Point>(*points), 3);
    std::vector<Point>* grahamOmpResult4Threads =
        grahamOmp(new std::vector<Point>(*points), 4);
    std::vector<Point>* grahamOmpResult8Threads =
        grahamOmp(new std::vector<Point>(*points), 8);
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult2Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult3Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult4Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult8Threads, *points));
}

TEST(Graham, Test_Random_1000_Points_In_Range_1000_2000) {
    auto points = generateRandomPointsVector(1000, 1000, 2000);
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    std::vector<Point>* grahamOmpResult2Threads =
        grahamOmp(new std::vector<Point>(*points), 2);
    std::vector<Point>* grahamOmpResult3Threads =
        grahamOmp(new std::vector<Point>(*points), 3);
    std::vector<Point>* grahamOmpResult4Threads =
        grahamOmp(new std::vector<Point>(*points), 4);
    std::vector<Point>* grahamOmpResult8Threads =
        grahamOmp(new std::vector<Point>(*points), 8);
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult2Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult3Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult4Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult8Threads, *points));
}

TEST(Graham, Test_Random_1_000_000_Points_In_Range_10000_20000) {
    auto points = generateRandomPointsVector(1000000, 10000, 20000);
    std::vector<Point>* grahamResult =
        graham(new std::vector<Point>(*points), 0, points->size());
    std::vector<Point>* grahamOmpResult2Threads =
        grahamOmp(new std::vector<Point>(*points), 2);
    std::vector<Point>* grahamOmpResult3Threads =
        grahamOmp(new std::vector<Point>(*points), 3);
    std::vector<Point>* grahamOmpResult4Threads =
        grahamOmp(new std::vector<Point>(*points), 4);
    std::vector<Point>* grahamOmpResult8Threads =
        grahamOmp(new std::vector<Point>(*points), 8);
    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult2Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult3Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult4Threads, *points));
    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult8Threads, *points));
}

// This test takes a lot of time to execute (~80 seconds),
// so should be used only for demonstration purposes

// TEST(Graham, Test_Random_10_000_000_Points_In_Range_100000_200000) {
//    auto points = generateRandomPointsVector(10000000, 100000, 200000);
//    std::vector<Point>* grahamResult =
//        graham(new std::vector<Point>(*points), 0, points->size());
//    std::vector<Point>* grahamOmpResult2Threads =
//        grahamOmp(new std::vector<Point>(*points), 2);
//    std::vector<Point>* grahamOmpResult3Threads =
//        grahamOmp(new std::vector<Point>(*points), 3);
//    std::vector<Point>* grahamOmpResult4Threads =
//        grahamOmp(new std::vector<Point>(*points), 4);
//    std::vector<Point>* grahamOmpResult8Threads =
//        grahamOmp(new std::vector<Point>(*points), 8);
//    ASSERT_TRUE(isGrahamCorrect(*grahamResult, *points));
//    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult2Threads, *points));
//    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult3Threads, *points));
//    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult4Threads, *points));
//    ASSERT_TRUE(isGrahamCorrect(*grahamOmpResult8Threads, *points));
//}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
