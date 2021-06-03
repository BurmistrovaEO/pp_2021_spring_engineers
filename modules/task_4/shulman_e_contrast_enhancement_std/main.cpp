// Copyright 2021 Shulman Egor
#include <gtest/gtest.h>
#include <time.h>
#include "./contrast_enhancement_std.h"

TEST(STD, testCreateMatrix) {
    int n = 100, m = 100;
    ASSERT_NO_THROW(getRandomMatrix(n, m));
}

TEST(STD, testRowsNegative) {
    int n = -100, m = 100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(TBB, testColumnsNegative) {
    int n = 100, m = -100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(STD, testRowsAndColumnsNegative) {
    int n = -100, m = -100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(STD, testContrastEnhancementSeq) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancementSeq(matrix));
}

TEST(STD, testContrastEnhancementStd) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancementStd(matrix));
}

TEST(STD, testComparisonSeqWithExpected) {
    int n = 5, m = 5;

    std::vector<int> matrix(n * m);
    std::vector<int> result(n * m);
    for (int i = 0; i < n * m; ++i) {
        matrix[i] = i + 10;
        result[i] = i + 10;
    }

    matrix = ContrastEnhancementSeq(matrix);

    int yMax = *std::max_element(result.begin(), result.end());
    int yMin = *std::min_element(result.begin(), result.end());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = ((result[i] - yMin) * 255) / (yMax - yMin);

    ASSERT_EQ(result, matrix);
}

TEST(STD, testComparisonStdWithExpected) {
    int n = 5, m = 5;

    std::vector<int> matrix(n * m);
    std::vector<int> result(n * m);
    for (int i = 0; i < n * m; ++i) {
        matrix[i] = i + 10;
        result[i] = i + 10;
    }

    matrix = ContrastEnhancementStd(matrix);

    int yMax = *std::max_element(result.begin(), result.end());
    int yMin = *std::min_element(result.begin(), result.end());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = ((result[i] - yMin) * 255) / (yMax - yMin);

    ASSERT_EQ(result, matrix);
}

TEST(STD, testComparisonStdWithSeq) {
    int n = 50, m = 50;
    std::vector<int> matrix = getRandomMatrix(n, m);
    std::vector<int> resultSeq(n * m);
    std::vector<int> resultStd(n * m);

    clock_t startSeq = clock();
    resultSeq = ContrastEnhancementSeq(matrix);
    clock_t endSeq = clock();
    double timeSeq = static_cast<double>(endSeq - startSeq) / CLOCKS_PER_SEC;

    clock_t startStd = clock();
    resultStd = ContrastEnhancementStd(matrix);
    clock_t endStd = clock();
    double timeStd = static_cast<double>(endStd - startStd) / CLOCKS_PER_SEC;

    printf("Sequential: %f seconds\n", timeSeq);
    printf("Std:        %f seconds\n", timeStd);

    ASSERT_EQ(resultSeq, resultStd);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
