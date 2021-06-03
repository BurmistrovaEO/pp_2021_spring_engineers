// Copyright 2021 Shulman Egor
#include <omp.h>
#include <gtest/gtest.h>
#include "./contrast_enhancement_omp.h"

TEST(OMP, testCreateMatrix) {
    int n = 100, m = 100;
    ASSERT_NO_THROW(getRandomMatrix(n, m));
}

TEST(OMP, testRowsNegative) {
    int n = -100, m = 100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(OMP, testColumnsNegative) {
    int n = 100, m = -100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(OMP, testRowsAndColumnsNegative) {
    int n = -100, m = -100;
    ASSERT_ANY_THROW(getRandomMatrix(n, m));
}

TEST(OMP, testContrastEnhancement) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancement(matrix));
}

TEST(OMP, testContrastEnhancementOmp) {
    int n = 100, m = 100;
    std::vector<int> matrix = getRandomMatrix(n, m);
    ASSERT_NO_THROW(ContrastEnhancementOMP(matrix));
}

TEST(OMP, testComparisonSeqWithExpected) {
    int n = 5, m = 5;

    std::vector<int> matrix(n * m);
    std::vector<int> result(n * m);
    for (int i = 0; i < n * m; ++i) {
        matrix[i] = i + 10;
        result[i] = i + 10;
    }

    matrix = ContrastEnhancement(matrix);

    int yMax = *std::max_element(result.begin(), result.end());
    int yMin = *std::min_element(result.begin(), result.end());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = (result[i] - yMin) * (255 / (yMax - yMin));

    ASSERT_EQ(result, matrix);;
}

TEST(OMP, testComparisonOmpWithExpected) {
    int n = 5, m = 5;

    std::vector<int> matrix(n * m);
    std::vector<int> result(n * m);
    for (int i = 0; i < n * m; ++i) {
        matrix[i] = i + 10;
        result[i] = i + 10;
    }

    matrix = ContrastEnhancementOMP(matrix);

    int yMax = *std::max_element(result.begin(), result.end());
    int yMin = *std::min_element(result.begin(), result.end());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = (result[i] - yMin) * (255 / (yMax - yMin));

    ASSERT_EQ(result, matrix);;
}

TEST(OMP, testComparisonSeqWithOmp) {
    int n = 5000, m = 5000;
    std::vector<int> matrix = getRandomMatrix(n, m);
    std::vector<int> resultSeq(n * m);
    std::vector<int> resultOmp(n * m);

    double startSeq = 0., endSeq = 0., startOmp = 0., endOmp = 0.;

    startSeq = omp_get_wtime();
    resultSeq = ContrastEnhancement(matrix);
    endSeq = omp_get_wtime();

    startOmp = omp_get_wtime();
    resultOmp = ContrastEnhancementOMP(matrix);
    endOmp = omp_get_wtime();

    double timeSeq = endSeq - startSeq;
    double timeOmp = endOmp - startOmp;
    double tick = omp_get_wtick();

    std::cout << "Sequential: " << timeSeq
              << " Omp: "       << timeOmp
              << " Precision: " << tick
              << std::endl;

    ASSERT_EQ(resultSeq, resultOmp);;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
