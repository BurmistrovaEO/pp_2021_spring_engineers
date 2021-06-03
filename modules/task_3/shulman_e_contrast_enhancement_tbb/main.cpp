// Copyright 2021 Shulman Egor
#include <omp.h>
#include <gtest/gtest.h>
#include <tbb/tick_count.h>
#include "./contrast_enhancement_tbb.h"

TEST(TBB, testCreateMatrix) {
    int n = 100, m = 100;
    ASSERT_NO_THROW(getRandomMatrix(n, m));
}

TEST(TBB, testRowsNegative) {
    int n = -100, m = 100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(TBB, testColumnsNegative) {
    int n = 100, m = -100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(TBB, testRowsAndColumnsNegative) {
    int n = -100, m = -100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(TBB, testContrastEnhancement) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancementSeq(matrix));
}

TEST(TBB, testContrastEnhancementOmp) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancementOmp(matrix));
}

TEST(TBB, testContrastEnhancementTbb) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancementTbb(matrix));
}

TEST(TBB, testComparisonSeqWithExpected) {
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

    ASSERT_EQ(result, matrix);;
}

TEST(TBB, testComparisonOmpWithExpected) {
    int n = 5, m = 5;

    std::vector<int> matrix(n * m);
    std::vector<int> result(n * m);
    for (int i = 0; i < n * m; ++i) {
        matrix[i] = i + 10;
        result[i] = i + 10;
    }

    matrix = ContrastEnhancementOmp(matrix);

    int yMax = *std::max_element(result.begin(), result.end());
    int yMin = *std::min_element(result.begin(), result.end());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = ((result[i] - yMin) * 255) / (yMax - yMin);

    ASSERT_EQ(result, matrix);;
}

TEST(TBB, testComparisonTbbWithExpected) {
    int n = 5, m = 5;

    std::vector<int> matrix(n * m);
    std::vector<int> result(n * m);
    for (int i = 0; i < n * m; ++i) {
        matrix[i] = i + 10;
        result[i] = i + 10;
    }

    matrix = ContrastEnhancementTbb(matrix);

    int yMax = *std::max_element(result.begin(), result.end());
    int yMin = *std::min_element(result.begin(), result.end());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = ((result[i] - yMin) * 255) / (yMax - yMin);

    ASSERT_EQ(result, matrix);;
}

TEST(TBB, testComparisonOmpWithTbb) {
    int n = 1000, m = 1000;
    std::vector<int> matrix = getRandomMatrix(n, m);
    std::vector<int> resultSeq(n * m);
    std::vector<int> resultOmp(n * m);
    std::vector<int> resultTbb(n * m);

    tbb::tick_count startSeq = tbb::tick_count::now();
    resultSeq = ContrastEnhancementSeq(matrix);
    tbb::tick_count endSeq = tbb::tick_count::now();

    tbb::tick_count startOmp = tbb::tick_count::now();
    resultOmp = ContrastEnhancementOmp(matrix);
    tbb::tick_count endOmp = tbb::tick_count::now();

    tbb::tick_count startTbb = tbb::tick_count::now();
    resultTbb = ContrastEnhancementTbb(matrix);
    tbb::tick_count endTbb = tbb::tick_count::now();

    printf("Sequential: %f seconds\n", (endSeq - startSeq).seconds());
    printf("Omp:        %f seconds\n", (endOmp - startOmp).seconds());
    printf("Tbb:        %f seconds\n", (endTbb - startTbb).seconds());

    ASSERT_EQ(resultOmp, resultTbb);;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
