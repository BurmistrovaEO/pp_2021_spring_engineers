// Copyright 2021 Ognev Denis

#include <../../modules/task_2/ognev_d_markings_hulls_omp/markings_hulls_omp.h>
#include <omp.h>

#include <algorithm>
#include <ctime>
#include <iostream>
#include <list>
#include <random>
#include <set>
#include <utility>
#include <vector>

int border = 135;

std::vector<int> treeUsual;
std::vector<int> treeOMP;
std::vector<int> treeOMP2;

std::vector<int> hull;
std::vector<int> hullOMP;

std::vector<bool> GenerateRandomMap(int width, int height) {
  if (width < 1 || height < 1) throw;
  std::mt19937 random;
  std::uniform_real_distribution<> urd(0, 1);
  random.seed(static_cast<unsigned int>(time(0)));
  std::vector<bool> result(width * height);

  for (int i = 0; i < width * height; i++) {
    float r = urd(random);
    if (r > 0.7f)
      result[i] = false;
    else
      result[i] = true;
  }

  return result;
}

std::vector<bool> GenerateConnectedMap(int width, int height) {
  if (width < 1 || height < 1) throw;
  std::mt19937 random;
  std::uniform_real_distribution<> urd(0, 1);
  random.seed(static_cast<unsigned int>(time(0)));
  std::vector<bool> result(width * height);
  int index;
  float r1 = urd(random);
  if (r1 > 0.5f)
    result[0] = true;
  else
    result[0] = false;
  for (int i = 1; i < width * height; i++) {
    float r = urd(random);
    if (i > width) {
      if (i % width != 0) {
        index = i - width;
        if (result[index] == true) {
          index = i - 1;
          if (result[index] == true) {
            if (r > 0.2f) {
              result[i] = true;
            } else {
              result[i] = false;
            }
          } else {
            if (r > 0.4f) {
              result[i] = true;
            } else {
              result[i] = false;
            }
          }
        } else {
          index = i - 1;
          if (result[index] == true) {
            if (r > 0.4f) {
              result[i] = false;
            } else {
              result[i] = true;
            }
          } else {
            if (r > 0.2f) {
              result[i] = false;
            } else {
              result[i] = true;
            }
          }
        }
      }
    } else {
      index = i - 1;
      if (result[index] == true) {
        if (r > 0.2f) {
          result[i] = true;
        } else {
          result[i] = false;
        }
      } else {
        if (r > 0.2f) {
          result[i] = false;
        } else {
          result[i] = true;
        }
      }
      continue;
    }
    if (i % width == 0) {
      index = i - width;
      if (result[index] == true) {
        if (r > 0.2f) {
          result[i] = true;
        } else {
          result[i] = false;
        }
      } else {
        if (r > 0.2f) {
          result[i] = false;
        } else {
          result[i] = true;
        }
      }
    }
  }

  return result;
}

int findInTree(int x) {
  while (treeUsual[x] != 0) x = treeUsual[x];
  return x;
}

int findInTreeOMP(int x) {
  while (treeOMP[x] != 0) x = treeOMP[x];
  return x;
}

int findInTreeOMPver2(int x) {
  while (treeOMP2[x] != 0) x = treeOMP2[x];
  return x;
}

void unionInTree(int x, int y) {
  x = findInTree(x);
  y = findInTree(y);
  if (y != x) treeUsual[y] = x;
}

void unionInTreeOMP(int x, int y) {
#pragma omp parallel
  {
#pragma omp sections nowait
    {
#pragma omp section
      { x = findInTreeOMP(x); }
#pragma omp section
      { y = findInTreeOMP(y); }
    }
  }
  if (y != x) treeOMP[y] = x;
}

void unionInTreeOMPver2(int x, int y) {
#pragma omp parallel
  {
#pragma omp sections nowait
    {
#pragma omp section
      { x = findInTreeOMPver2(x); }
#pragma omp section
      { y = findInTreeOMPver2(y); }
    }
  }
  if (y != x) treeOMP2[y] = x;
}

std::vector<int> firstStep(
    std::vector<bool> map, int width,
    int height) {  // use tree table for solving border problem
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  std::vector<int> result(width * height);
  int label = 0;
  treeUsual.clear();
  treeUsual = std::vector<int>(1, 0);

  for (int i = 0; i < width * height; i++) {
    result[i] = 0;
    if (map[i] != false) {
      if (label == 0) {
        treeUsual.push_back(0);
        label++;
        result[i] = label;
      } else if (i < width) {
        if (result[i - 1] != 0) {
          result[i] = result[i - 1];
        } else {
          treeUsual.push_back(0);
          label++;
          result[i] = label;
        }
      } else if (i % width == 0) {
        if (result[i - width] != 0) {
          result[i] = result[i - width];
        } else {
          treeUsual.push_back(0);
          label++;
          result[i] = label;
        }
      } else {
        if (result[i - width] != 0) {
          if (result[i - 1] != 0) {
            if (result[i - width] != result[i - 1]) {
              int first = result[i - width];
              int second = result[i - 1];
              unionInTree(first, second);
            }
          }
          result[i] = result[i - width];
        } else if (result[i - 1] != 0) {
          result[i] = result[i - 1];
        } else {
          treeUsual.push_back(0);
          label++;
          result[i] = label;
        }
      }
    }
  }
  return result;
}

std::vector<int> firstStepOMP(
    std::vector<bool> map, int width,
    int height) {  // use tree table for solving border problem
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  std::vector<int> result(width * height, 0);
  int label = 0;
  treeOMP.clear();
  treeOMP = std::vector<int>(1, 0);

  if (map[0] == true) {
    treeOMP.push_back(0);
    label++;
    result[0] = label;
  }

  int smallerDimsension = width > height ? height : width;

  for (int i = 0; i < smallerDimsension; i++) {
#pragma omp parallel
    {
#pragma omp single nowait
      for (int j = i; j < width; j++) {
        if (map[i * width + j] == true) {
          if (i == 0) {
            if (j == 0) continue;
            if (map[i * width + j - 1] == true) {
              result[i * width + j] = result[i * width + j - 1];
            } else {
#pragma omp critical
              treeOMP.push_back(0);
#pragma omp atomic
              label++;
              result[i * width + j] = label;
            }
          } else if (map[(i - 1) * width + j] == true) {
            if (map[i * width + j - 1] == true) {
              if (result[(i - 1) * width + j] != result[i * width + j - 1]) {
                int first = result[(i - 1) * width + j];
                int second = result[i * width + j - 1];
#pragma omp critical
                unionInTreeOMP(first, second);
              }
              result[i * width + j] = result[i * width + j - 1];
            } else {
              result[i * width + j] = result[(i - 1) * width + j];
            }
          } else if (map[i * width + j - 1] == true) {
            result[i * width + j] = result[i * width + j - 1];
          } else {
#pragma omp critical
            treeOMP.push_back(0);
#pragma omp atomic
            label++;
            result[i * width + j] = label;
          }
        }
      }

#pragma omp single nowait
      for (int j = i; j < height; j++) {
        if (map[j * width + i] == true) {
          if (i == 0) {
            if (j == 0) continue;
            if (map[(j - 1) * width + i] == true) {
              result[j * width + i] = result[(j - 1) * width + i];
            } else {
#pragma omp critical
              treeOMP.push_back(0);
#pragma omp critical
              label++;
              result[j * width + i] = label;
            }
          } else if (map[(j - 1) * width + i] == true) {
            if (map[j * width + i - 1] == true) {
              if (result[(j - 1) * width + i] != result[j * width + i - 1]) {
                int first = result[(j - 1) * width + i];
                int second = result[j * width + i - 1];
#pragma omp critical
                unionInTreeOMP(first, second);
              }
              result[j * width + i] = result[j * width + i - 1];
            } else {
              result[j * width + i] = result[(j - 1) * width + i];
            }
          } else if (map[j * width + i - 1] == true) {
            result[j * width + i] = result[j * width + i - 1];
          } else {
#pragma omp critical
            treeOMP.push_back(0);
#pragma omp atomic
            label++;
            result[j * width + i] = label;
          }
        }
      }
    }
  }
  return result;
}
std::vector<int> firstStepOMPver2(
    std::vector<bool> map, int width,
    int height) {  // use tree table for solving border problem
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  std::vector<int> result(width * height, 0);
  int label = 0;
  treeOMP2.clear();
  treeOMP2 = std::vector<int>(1, 0);

  for (int i = 0; i < width; i++) {
#pragma omp parallel for
    for (int j = 0; j <= i; j++) {
      if (j >= height) {
        continue;
      }

      if (map[j * width + (i - j)] == true) {
        if (j == 0) {
          if (i == 0) {
#pragma omp critical
            treeOMP2.push_back(0);
#pragma omp atomic
            label++;
            result[j * width + (i - j)] = label;
          } else if (map[j * width + (i - j) - 1] == true) {
            result[j * width + (i - j)] = result[j * width + (i - j) - 1];
          } else {
#pragma omp critical
            treeOMP2.push_back(0);
#pragma omp atomic
            label++;
            result[j * width + (i - j)] = label;
          }
        } else if (j == i) {
          if (map[j * width + (i - j) - width] == true) {
            result[j * width + (i - j)] = result[j * width + (i - j) - width];
          } else {
#pragma omp critical
            treeOMP2.push_back(0);
#pragma omp atomic
            label++;
            result[j * width + (i - j)] = label;
          }
        } else if (map[j * width + (i - j) - 1] == true) {
          if (map[j * width + (i - j) - width] == true) {
            if (result[j * width + (i - j) - width] !=
                result[j * width + (i - j) - 1]) {
              int first = result[j * width + (i - j) - width];
              int second = result[j * width + (i - j) - 1];
#pragma omp critical
              unionInTreeOMPver2(first, second);
            }
          }
          result[j * width + (i - j)] = result[j * width + (i - j) - 1];
        } else if (map[j * width + (i - j) - width] == true) {
          result[j * width + (i - j)] = result[j * width + (i - j) - width];
        } else {
#pragma omp critical
          treeOMP2.push_back(0);
#pragma omp atomic
          label++;
          result[j * width + (i - j)] = label;
        }
      }
    }
  }

  for (int i = width * 2 - 1; i < width * height; i += width) {
#pragma omp parallel for
    for (int j = 0; j < width * height - i; j += width) {
      if (j / width >= width) {
        continue;
      }

      if (map[(i + j) - j / width] == true) {
        if (j == width * (width - 1)) {
          if (map[(i + j) - j / width - width] == true) {
            result[(i + j) - j / width] = result[(i + j) - j / width - width];
          } else {
#pragma omp critical
            treeOMP2.push_back(0);
#pragma omp atomic
            label++;
            result[(i + j) - j / width] = label;
          }
        } else if (map[(i + j) - j / width - 1] == true) {
          if (map[(i + j) - j / width - width] == true) {
            if (result[(i + j) - j / width - width] !=
                result[(i + j) - j / width - 1]) {
              int first = result[(i + j) - j / width - 1];
              int second = result[(i + j) - j / width - width];
#pragma omp critical
              unionInTreeOMPver2(first, second);
            }
          }
          result[(i + j) - j / width] = result[(i + j) - j / width - 1];
        } else if (map[(i + j) - j / width - width] == true) {
          result[(i + j) - j / width] = result[(i + j) - j / width - width];
        } else {
#pragma omp critical
          treeOMP2.push_back(0);
#pragma omp atomic
          label++;
          result[(i + j) - j / width] = label;
        }
      }
    }
  }

  return result;
}

std::vector<int> secondStep(std::vector<int> map, int width, int height) {
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  if (treeUsual.size() ==
      1)  // if no elements were added to tree table -> no collisions
    return map;
  std::vector<int> result(map);
  for (int i = 0; i < width * height; i++) {
    if (result[i] != 0) {
      result[i] = findInTree(result[i]);  // searching parent
    }
  }
  return result;
}

std::vector<int> secondStepOMP(std::vector<int> map, int width, int height) {
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  if (treeOMP.size() ==
      1)  // if no elements were added to tree table -> no collisions
    return map;
  std::vector<int> result(map);
#pragma omp parallel for
  for (int i = 0; i < width * height; i++) {
    if (result[i] != 0) {
      result[i] = findInTreeOMP(result[i]);  // searching parent
    }
  }
  return result;
}

std::vector<int> secondStepOMPver2(std::vector<int> map, int width,
                                   int height) {
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  if (treeOMP2.size() ==
      1)  // if no elements were added to tree table -> no collisions
    return map;
  std::vector<int> result(map);
#pragma omp parallel for
  for (int i = 0; i < width * height; i++) {
    if (result[i] != 0) {
      result[i] = findInTreeOMPver2(result[i]);  // searching parent
    }
  }
  return result;
}

int orientTriangle2(int x1, int y1, int x2, int y2, int x3, int y3) {
  return (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2);
}

int orientTriangle(int x1, int y1, int x2, int y2, int x3, int y3) {
  return (x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1);
}

bool isInside(int x1, int y1, int x2, int y2, int x3, int y3) {
  if ((x2 >= x1 && x2 <= x3 && y2 >= y1 && y2 <= y3) ||
      (x2 >= x3 && x2 <= x1 && y2 >= y1 && y2 <= y3) ||
      (x2 >= x1 && x2 <= x3 && y2 >= y3 && y2 <= y1) ||
      (x2 >= x3 && x2 <= x1 && y2 >= y3 && y2 <= y1)) {
    return true;
  }
  return false;
}

std::vector<int> jarvis(std::vector<int> map, int width, int height, int mark,
                        int index) {
  if (width < 1 || height < 1 || index < 0 || mark < 0 || map.empty()) {
    throw "";
  }
  std::list<int> vertecies;
  std::vector<int> temp(map);
  int first = index;
  int current = index;

  map[index] = border;

  do {
    bool marked = false;
    int next = current;
    for (int i = 0; i < width * height; i++) {
      if (temp[i] == mark) {
        int sign =
            orientTriangle2(current / width, current % width, next / width,
                            next % width, i / width, i % width);
        if (sign > 0) {
          marked = true;
          next = i;
        } else if (sign == 0) {
          if (isInside(current / width, current % width, next / width,
                       next % width, i / width, i % width)) {
            marked = true;
            next = i;
          }
        }
      }
    }
    if (marked && next != first) {
      marked = false;
      map[next] = border;
    }
    current = next;
  } while (current != first);

  border++;
  return map;
}

std::vector<int> convexHull(std::vector<int> map, int width, int height) {
  if (width < 1 || height < 1 || map.empty()) {
    throw "";
  }
  std::vector<int> result(map);
  std::set<int> marks(map.begin(), map.end());
  marks.erase(0);

  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      if (result[i * width + j] != 0) {
        if (marks.find(result[i * width + j]) != marks.end()) {
          marks.erase(result[i * width + j]);
          result = jarvis(result, width, height, result[i * width + j],
                          i * width + j);
        }
      }
    }
  }

  return result;
}

void findHull(std::vector<int> v, int leftX, int leftY, int rightX, int rightY,
              int maxIndexX, int maxIndexY, int width, int height) {
  if (maxIndexY * width + maxIndexX >= 0) {
    hull.push_back(maxIndexY * width + maxIndexX);
  }
  if (v.size() == 0) return;

  std::vector<int> S1;
  std::vector<int> S2;
  int leftMaxIndexX = -1, rightMaxIndexX = -1;
  int leftMaxIndexY = -1, rightMaxIndexY = -1;
  int leftMax = 0, rightMax = 0;

  int size = v.size();
  for (int i = 0; i < size; i++) {
    int signL = orientTriangle(leftX, leftY, maxIndexX, maxIndexY, v[i] % width,
                               v[i] / width);
    int signR = orientTriangle(maxIndexX, maxIndexY, rightX, rightY,
                               v[i] % width, v[i] / width);
    if (signL > 0) {
      if (signL > leftMax) {
        leftMax = signL;
        leftMaxIndexX = v[i] % width;
        leftMaxIndexY = v[i] / width;
      }
      S1.push_back(v[i]);
    }
    if (signR > 0) {
      if (signR > rightMax) {
        rightMax = signR;
        rightMaxIndexX = v[i] % width;
        rightMaxIndexY = v[i] / width;
      }
      S2.push_back(v[i]);
    }
  }

  {
    if (leftMax != 0)
      findHull(S1, leftX, leftY, maxIndexX, maxIndexY, leftMaxIndexX,
               leftMaxIndexY, width, height);
    if (rightMax != 0)
      findHull(S2, maxIndexX, maxIndexY, rightX, rightY, rightMaxIndexX,
               rightMaxIndexY, width, height);
  }
}

std::vector<int> quickHull(std::vector<int> map, int width, int height,
                           int mark, int leftPos, int rightPos) {
  int leftX = leftPos % width;
  int leftY = leftPos / width;
  int rightX = rightPos % width;
  int rightY = rightPos / width;
  std::vector<int> result(map);
  hull.clear();
  hull.push_back(leftPos);
  hull.push_back(rightPos);
  std::vector<int> S1;
  std::vector<int> S2;
  int leftMaxIndexX = -1, rightMaxIndexX = -1;
  int leftMaxIndexY = -1, rightMaxIndexY = -1;
  int leftMax = 0, rightMax = 0;

  for (int i = leftPos; i < rightPos; i++) {
    if (map[i] == mark) {
      int sign =
          orientTriangle(leftX, leftY, rightX, rightY, i % width, i / width);
      if (sign > 0) {
        if (sign > leftMax) {
          leftMax = sign;
          leftMaxIndexX = i % width;
          leftMaxIndexY = i / width;
        }
        S1.push_back(i);
      } else if (sign < 0) {
        if (sign < rightMax) {
          rightMax = sign;
          rightMaxIndexX = i % width;
          rightMaxIndexY = i / width;
        }
        S2.push_back(i);
      }
    }
  }

  if (leftMax != 0)
    findHull(S1, leftX, leftY, rightX, rightY, leftMaxIndexX, leftMaxIndexY,
             width, height);

  if (rightMax != 0)
    findHull(S2, rightX, rightY, leftX, leftY, rightMaxIndexX, rightMaxIndexY,
             width, height);

  int convSize = hull.size();
  for (int i = 0; i < convSize; i++) {
    result[hull[i]] = border;
  }
  border++;

  return result;
}

std::vector<int> convexQuickHull(std::vector<int> map, int width, int height) {
  std::vector<int> result(map);
  std::set<int> marks(map.begin(), map.end());
  marks.erase(0);

  for (int i = 0; i < width * height; i++) {
    int mark = map[i];
    if (mark != 0) {
      if (marks.find(mark) != marks.end()) {
        marks.erase(mark);
        for (int j = width * height - 1; j >= i; j--) {
          if (map[j] == mark) {
            result = quickHull(result, width, height, mark, i, j);
            break;
          }
        }
      }
    }
  }

  return result;
}

void findHullOMP(std::vector<int> v, std::vector<int>* hull, int leftX,
                 int leftY, int rightX, int rightY, int maxIndexX,
                 int maxIndexY, int width, int height) {
  // #pragma omp critical
  //  {
  if (maxIndexY * width + maxIndexX >= 0)
    (*hull).push_back(maxIndexY * width + maxIndexX);
  //    }
  //  }
  if (v.size() == 0) return;
  int size = v.size();

  std::vector<int> S1;
  std::vector<int> S2;

  int leftMaxIndexX = -1, rightMaxIndexX = -1;
  int leftMaxIndexY = -1, rightMaxIndexY = -1;
  int leftMax = 0, rightMax = 0;

  for (int i = 0; i < size; i++) {
    int signL = orientTriangle(leftX, leftY, maxIndexX, maxIndexY, v[i] % width,
                               v[i] / width);
    int signR = orientTriangle(maxIndexX, maxIndexY, rightX, rightY,
                               v[i] % width, v[i] / width);
    if (signL > 0) {
      if (signL > leftMax) {
        leftMax = signL;
        leftMaxIndexX = v[i] % width;
        leftMaxIndexY = v[i] / width;
      }
      S1.push_back(v[i]);
    }
    if (signR > 0) {
      if (signR > rightMax) {
        rightMax = signR;
        rightMaxIndexX = v[i] % width;
        rightMaxIndexY = v[i] / width;
      }
      S2.push_back(v[i]);
    }
  }

  // #pragma omp parallel
  {
    // #pragma omp sections nowait
    {
      // #pragma omp section
      if (leftMax != 0)
        findHullOMP(S1, hull, leftX, leftY, maxIndexX, maxIndexY, leftMaxIndexX,
                    leftMaxIndexY, width, height);
      // #pragma omp section
      if (rightMax != 0)
        findHullOMP(S2, hull, maxIndexX, maxIndexY, rightX, rightY,
                    rightMaxIndexX, rightMaxIndexY, width, height);
    }
  }
}

std::vector<int> quickHullOMP(std::vector<int> map, int width, int height,
                              int mark, int leftPos, int rightPos) {
  int leftX = leftPos % width;
  int leftY = leftPos / width;
  int rightX = rightPos % width;
  int rightY = rightPos / width;

  std::vector<int> hull;

  hull.push_back(leftPos);
  hull.push_back(rightPos);

  std::vector<int> S1;
  std::vector<int> S2;

  int leftMaxIndexX = -1, rightMaxIndexX = -1;
  int leftMaxIndexY = -1, rightMaxIndexY = -1;
  int leftMax = 0, rightMax = 0;

  for (int i = leftPos; i < rightPos; i++) {
    if (map[i] == mark) {
      int sign =
          orientTriangle(leftX, leftY, rightX, rightY, i % width, i / width);
      // left
      if (sign > 0) {
        if (sign > leftMax) {
          leftMax = sign;
          leftMaxIndexX = i % width;
          leftMaxIndexY = i / width;
        }
        S1.push_back(i);
        // right
      } else if (sign < 0) {
        if (sign < rightMax) {
          rightMax = sign;
          rightMaxIndexX = i % width;
          rightMaxIndexY = i / width;
        }
        S2.push_back(i);
      }
    }
  }

  // #pragma omp parallel
  {
    // #pragma omp sections nowait
    {
      // #pragma omp section
      if (leftMax != 0)
        findHullOMP(S1, &hull, leftX, leftY, rightX, rightY, leftMaxIndexX,
                    leftMaxIndexY, width, height);

      // #pragma omp section
      if (rightMax != 0)
        findHullOMP(S2, &hull, rightX, rightY, leftX, leftY, rightMaxIndexX,
                    rightMaxIndexY, width, height);
    }
  }

  std::vector<int> result(hull);

  return result;
}

std::vector<int> convexQuickHullOMP(std::vector<int> map, int width,
                                    int height) {
  std::vector<int> result(map);
  std::set<int> marks(map.begin(), map.end());
  marks.erase(0);
  std::vector<std::pair<int, int>> vert;
  std::vector<int> _marks;

  for (int i = 0; i < width * height; i++) {
    int mark = map[i];
    if (mark != 0) {
      if (marks.find(mark) != marks.end()) {
        marks.erase(mark);
        for (int j = width * height - 1; j >= i; j--) {
          if (map[j] == mark) {
            vert.push_back(std::pair<int, int>(i, j));
            _marks.push_back(mark);
            break;
          }
        }
      }
    }
  }

  std::vector<std::vector<int>> hulls(vert.size());
  std::vector<int> borders(vert.size());

  int vertSize = vert.size();
  int hullsSize = hulls.size();

#pragma omp parallel for
  for (int i = 0; i < vertSize; i++) {
    hulls[i] = quickHullOMP(map, width, height, _marks[i], vert[i].first,
                            vert[i].second);
  }

  for (int i = 0; i < hullsSize; i++) {
    borders[i] = border++;
  }

#pragma omp parallel for
  for (int i = 0; i < hullsSize; i++) {
    int hullSize = hulls[i].size();
    for (int j = 0; j < hullSize; j++) {
      result[hulls[i][j]] = borders[i];
    }
  }

  return result;
}
