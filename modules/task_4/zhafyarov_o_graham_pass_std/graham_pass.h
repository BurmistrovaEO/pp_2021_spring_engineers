// Copyright 2021 Zhafyarov Oleg
#ifndef MODULES_TASK_4_ZHAFYAROV_O_GRAHAM_PASS_STD_GRAHAM_PASS_H_
#define MODULES_TASK_4_ZHAFYAROV_O_GRAHAM_PASS_STD_GRAHAM_PASS_H_

// #include <tbb/tbb.h>
#include <omp.h>
#include <vector>
#include <utility>

struct point{
  int x;
  int y;
};

bool CompareVectors(const std::vector<point>& vec1,
                    const std::vector<point>& vec2);
std::vector<point> RandomVector(int size);
double Rotation(point a, point b, point c);
std::vector<point> GrahamPassSeq(const std::vector<point>& basisVec);
/*
std::vector<point> GrahamPassOmp(const std::vector<point>& basisVec);
*/
/*
std::vector<point> GrahamPassTbb(
    const std::vector<point>& basisVec, int grainsize);
*/
void GrahamPassSeqForStd(const std::vector<point>& basis_vec,
                         std::vector<point>* result_after);
std::vector<point> GrahamPassStd(const std::vector<point>& basis_vec);
#endif  // MODULES_TASK_4_ZHAFYAROV_O_GRAHAM_PASS_STD_GRAHAM_PASS_H_
