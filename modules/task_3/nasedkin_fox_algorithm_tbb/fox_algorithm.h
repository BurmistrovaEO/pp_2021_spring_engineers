// Copyright 2021 Nasedkin Nikita
#ifndef MODULES_TASK_3_NASEDKIN_FOX_ALGORITHM_TBB_FOX_ALGORITHM_H_
#define MODULES_TASK_3_NASEDKIN_FOX_ALGORITHM_TBB_FOX_ALGORITHM_H_

#include <random>
#include <ctime>
#include <vector>

std::vector<double> MatrixInit(int size);
std::vector<double> SeqMatrixCalc(std::vector<double> a,
  std::vector<double> b, int size);
std::vector<double> SeqFox(std::vector<double> a,
  std::vector<double> b, int blockSize);
std::vector<double> TbbFox(const std::vector<double>& a,
  const std::vector<double>& b, int matrixSize, int countOfThreads);

#endif  // MODULES_TASK_3_NASEDKIN_FOX_ALGORITHM_TBB_FOX_ALGORITHM_H_
