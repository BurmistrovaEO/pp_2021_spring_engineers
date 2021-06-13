// Copyright 2021 Ekaterina Burmistrova
#include <omp.h>
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <functional>
#include "./trapezoid_integral_omp.h"

double integral(double x, double y = 1, double z = 1) {
    return (pow(x, 2));
}

double integral2(double x, double y, double z = 1) {
    return y / x;
}

double integral3(double x, double y, double z) {
    return x * y * z;
}

double integral4(double x, double y, double z = 1) {
    return pow(x, 2) + pow(y, 3);
}

double integral5(double x, double y, double z) {
    return x + y + z;
}

TEST(Parallel_Operations_OpenMP, Test_OneDim) {
    std::vector<std::pair<int, int>> vec(1);
    vec[0] = std::make_pair(2, 5);
    // std::function<double(double)> func= pow(x, 2);
    std::function<double(double, double, double)> ptr1 = integral;
    double par_int = SolveParallel(vec, ptr1);
    // int parallel_sum = getParallelOperations(vec, "+");
    vec.clear();
    ASSERT_NEAR(39.0, par_int, 0.6);
}

TEST(Parallel_Operations_OpenMP, Test_TwoDim) {
    std::vector<std::pair<int, int>> vec(2);
    vec[0] = std::make_pair(1, 2);
    vec[1] = std::make_pair(2, 3);
    // std::function<double(double)> func= pow(x, 2);
    double t1 = 0, t2 = 0;
    std::function<double(double, double, double)> ptr2 = integral2;
    t1 = omp_get_wtime();
    double par_int = SolveParallel(vec, ptr2);
    t2 = omp_get_wtime();
    std::cout << "Test runtime = " << t2 - t1 <<"\n";
    // int parallel_sum = getParallelOperations(vec, "+");
    vec.clear();
    ASSERT_NEAR(1.73, par_int, 0.5);
}

TEST(Parallel_Operations_OpenMP, Test_ThreeDim) {
    std::vector<std::pair<int, int>> vec(3);
    vec[0] = std::make_pair(1, 2);
    vec[1] = std::make_pair(2, 3);
    vec[2] = std::make_pair(4, 5);
    // std::function<double(double)> func= pow(x, 2);
    double t1 = 0, t2 = 0;
    std::function<double(double, double, double)> ptr3 = integral3;
    t1 = omp_get_wtime();
    double par_int = SolveParallel(vec, ptr3);
    t2 = omp_get_wtime();
    std::cout << "Test runtime = " << t2 - t1 <<"\n";
    // int parallel_sum = getParallelOperations(vec, "+");
    vec.clear();
    ASSERT_NEAR(16.87, par_int, 7.0);
}

TEST(Parallel_Operations_OpenMP, Test_TwoDimSum) {
    std::vector<std::pair<int, int>> vec(2);
    vec[0] = std::make_pair(1, 3);
    vec[1] = std::make_pair(2, 4);
    // std::function<double(double)> func= pow(x, 2);
    double t1 = 0, t2 = 0;
    std::function<double(double, double, double)> ptr4 = integral4;
    t1 = omp_get_wtime();
    double par_int = SolveParallelSum(vec, ptr4);
    t2 = omp_get_wtime();
    std::cout << "Test runtime = " << t2 - t1 <<"\n";
    // int parallel_sum = getParallelOperations(vec, "+");
    vec.clear();
    ASSERT_NEAR(137.33, par_int, 12.00);
}

TEST(Parallel_Operations_OpenMP, Test_ThreeDimSum) {
    std::vector<std::pair<int, int>> vec(3);
    vec[0] = std::make_pair(1, 4);
    vec[1] = std::make_pair(1, 2);
    vec[2] = std::make_pair(1, 3);
    // std::function<double(double)> func= pow(x, 2);
    double t1 = 0, t2 = 0;
    std::function<double(double, double, double)> ptr5 = integral5;
    t1 = omp_get_wtime();
    double par_int = SolveParallelSum(vec, ptr5);
    t2 = omp_get_wtime();
    std::cout<< "Test runtime = " << t2 - t1;
    // int parallel_sum = getParallelOperations(vec, "+");
    vec.clear();
    ASSERT_NEAR(36.00, par_int, 12.00);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
