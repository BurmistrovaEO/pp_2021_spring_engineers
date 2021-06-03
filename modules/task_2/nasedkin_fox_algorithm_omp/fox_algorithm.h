// Copyright 2021 Nasedkin Nikita
#ifndef MODULES_TASK_2_NASEDKIN_FOX_ALGORITHM_OMP_FOX_ALGORITHM_H_
#define MODULES_TASK_2_NASEDKIN_FOX_ALGORITHM_OMP_FOX_ALGORITHM_H_

#include <omp.h>
#include <random>
#include <ctime>
#include <vector>

std::vector<double> MatrixInit(int size);
std::vector<double> SeqMatrixCalc(std::vector<double> a,
  std::vector<double> b, int size);
std::vector<double> SeqFox(std::vector<double> a,
  std::vector<double> b, int blockSize);
std::vector<double> OmpFox(const std::vector<double>& a,
  const std::vector<double>& b, int matrixSize);

#endif  // MODULES_TASK_2_NASEDKIN_FOX_ALGORITHM_OMP_FOX_ALGORITHM_H_
