// Copyright 2021 Evseev Alexander
#include <gtest/gtest.h>
#include <omp.h>
#include <vector>
#include<iostream>
#include "./evseev_a_gauss_omp.h"

TEST(GaussFilter_Test, invalid_matrix) {
  int hight = 0;
  int widht = -1;
  ASSERT_ANY_THROW(RandMatrix(hight, widht));
}
TEST(OMP_GaussFilter_Test, Test_1500x1000x60thr) {
  int rows = 1500;
  int cols = 1000;
  double sigma = 7;
  int radius = 2;
  Matrix SEQ(rows, std::vector<double>(cols));
  Matrix OMP(rows, std::vector<double>(cols));
  Matrix mat(rows, std::vector<double>(cols));
  Matrix kernel = GaussKernel(radius, sigma);
  mat = RandMatrix(rows, cols);
  double start = omp_get_wtime();
  SEQ = SeqGauss(mat, rows, cols, radius, sigma);
  double end = omp_get_wtime();
  std::cout << "Sequential time: " << end - start << std::endl;
  start = omp_get_wtime();
  omp_set_num_threads(60);
  OMP = OmpGauss(mat, rows, cols, radius, sigma);
  end = omp_get_wtime();
  std::cout << "Parallel time: " << end - start << std::endl;
  ASSERT_EQ(SEQ, OMP);
}

TEST(OMP_GaussFilter_Test, Test_400x400_20thr) {
  int rows = 400;
  int cols = 400;
  double sigma = 7;
  int radius = 2;
  Matrix SEQ(rows, std::vector<double>(cols));
  Matrix OMP(rows, std::vector<double>(cols));
  Matrix mat(rows, std::vector<double>(cols));
  Matrix kernel = GaussKernel(radius, sigma);
  mat = RandMatrix(rows, cols);
  double start = omp_get_wtime();
  SEQ = SeqGauss(mat, rows, cols, radius, sigma);
  double end = omp_get_wtime();
  std::cout << "Sequential time: " << end - start << std::endl;
  start = omp_get_wtime();
  omp_set_num_threads(20);
  OMP = OmpGauss(mat, rows, cols, radius, sigma);
  end = omp_get_wtime();
  std::cout << "Parallel time: " << end - start << std::endl;
  ASSERT_EQ(SEQ, OMP);
}

TEST(GaussFilter_Test, Test_4x3) {
  int height = 4;
  int width = 3;
  double sigma = 1;
  int radius = 1;
  Matrix mat(height, std::vector<double>(width));
  Matrix gaussSeq(height, std::vector<double>(width));
  Matrix gaussOmp(height, std::vector<double>(width));
  Matrix resMat(height, std::vector<double>(width));
  Matrix KernelMat = GaussKernel(radius, sigma);
  mat[0][0] = 105;
  mat[0][1] = 177;
  mat[0][2] = 111;
  mat[1][0] = 200;
  mat[1][1] = 28;
  mat[1][2] = 151;
  mat[2][0] = 230;
  mat[2][1] = 44;
  mat[2][2] = 88;
  mat[3][0] = 112;
  mat[3][1] = 128;
  mat[3][2] = 155;

  resMat[0][0] = 98;
  resMat[0][1] = 170;
  resMat[0][2] = 104;
  resMat[1][0] = 194;
  resMat[1][1] = 25;
  resMat[1][2] = 146;
  resMat[2][0] = 228;
  resMat[2][1] = 38;
  resMat[2][2] = 81;
  resMat[3][0] = 109;
  resMat[3][1] = 122;
  resMat[3][2] = 147;

  gaussSeq = SeqGauss(mat, height, width, radius, sigma);
  gaussOmp = OmpGauss(mat, height, width, radius, sigma);
  ASSERT_EQ(resMat, gaussSeq);
  ASSERT_EQ(resMat, gaussOmp);
}

TEST(GaussFilter_Test, Test_8x1) {
  int height = 8;
  int width = 1;
  double sigma = 1;
  int radius = 1;
  Matrix mat(height, std::vector<double>(width));
  Matrix gaussSeq(height, std::vector<double>(width));
  Matrix gaussOmp(height, std::vector<double>(width));
  Matrix resMat(height, std::vector<double>(width));
  Matrix KernelMat = GaussKernel(radius, sigma);
  mat[0][0] = 21;
  mat[1][0] = 67;
  mat[2][0] = 2;
  mat[3][0] = 2;
  mat[4][0] = 77;
  mat[5][0] = 55;
  mat[6][0] = 12;
  mat[7][0] = 56;

  resMat[0][0] = 14;
  resMat[1][0] = 66;
  resMat[2][0] = 0;
  resMat[3][0] = 0;
  resMat[4][0] = 73;
  resMat[5][0] = 50;
  resMat[6][0] = 7;
  resMat[7][0] = 50;

  gaussSeq = SeqGauss(mat, height, width, radius, sigma);
  gaussOmp = OmpGauss(mat, height, width, radius, sigma);
  ASSERT_EQ(resMat, gaussSeq);
  ASSERT_EQ(resMat, gaussOmp);
}
