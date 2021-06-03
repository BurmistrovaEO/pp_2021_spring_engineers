// Copyright 2021 Ognev Denis
#ifndef MODULES_TASK_2_OGNEV_D_MARKINGS_HULLS_OMP_MARKINGS_HULLS_OMP_H_
#define MODULES_TASK_2_OGNEV_D_MARKINGS_HULLS_OMP_MARKINGS_HULLS_OMP_H_

#include <vector>

std::vector<bool> GenerateRandomMap(int width, int height);
std::vector<bool> GenerateConnectedMap(int width, int height);
int findInTree(int x);
int findInTreeOMP(int x);
int findInTreeOMPver2(int x);
void unionInTree(int x, int y);
void unionInTreeOMP(int x, int y);
void unionInTreeOMPver2(int x, int y);
std::vector<int> firstStep(std::vector<bool> map, int width, int height);
std::vector<int> firstStepOMP(std::vector<bool> map, int width, int height);
std::vector<int> firstStepOMPver2(std::vector<bool> map, int width, int height);
std::vector<int> secondStep(std::vector<int> map, int width, int height);
std::vector<int> secondStepOMP(std::vector<int> map, int width, int height);
std::vector<int> secondStepOMPver2(std::vector<int> map, int width, int height);
int orientTriangle2(int x1, int y1, int x2, int y2, int x3, int y3);
int orientTriangle(int x1, int y1, int x2, int y2, int x3, int y3);
bool isInside(int x1, int y1, int x2, int y2, int x3, int y3);
std::vector<int> jarvis(std::vector<int> map, int width, int height, int mark,
                        int index);
std::vector<int> convexHull(std::vector<int> map, int width, int height);
void findHull(std::vector<int> v, int leftX, int leftY, int rightX, int rightY,
              int maxIndexX, int maxIndexY, int width, int height);
std::vector<int> quickHull(std::vector<int> map, int width, int height,
                           int mark, int leftPos, int rightPos);
std::vector<int> convexQuickHull(std::vector<int> map, int width, int height);
void findHullOMP(std::vector<int> v, std::vector<int>* hull, int leftX,
                 int leftY, int rightX, int rightY, int maxIndexX,
                 int maxIndexY, int width, int height);
std::vector<int> quickHullOMP(std::vector<int> map, int width, int height,
                              int mark, int leftPos, int rightPos);
std::vector<int> convexQuickHullOMP(std::vector<int> map, int width,
                                    int height);

#endif  // MODULES_TASK_2_OGNEV_D_MARKINGS_HULLS_OMP_MARKINGS_HULLS_OMP_H_
