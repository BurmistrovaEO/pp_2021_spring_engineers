// Copyright 2021 Ognev Denis

#include <../../modules/task_3/ognev_d_markings_hulls_tbb/markings_hulls_tbb.h>
#include <tbb/tbb.h>

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

std::vector<int> hull;

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

void unionInTree(int x, int y) {
  x = findInTree(x);
  y = findInTree(y);
  if (y != x) treeUsual[y] = x;
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

void findHullTBB(std::vector<int> v, std::vector<int>* hull, int leftX,
                 int leftY, int rightX, int rightY, int maxIndexX,
                 int maxIndexY, int width, int height) {
  if (maxIndexY * width + maxIndexX >= 0)
    (*hull).push_back(maxIndexY * width + maxIndexX);

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

  if (leftMax != 0)
    findHullTBB(S1, hull, leftX, leftY, maxIndexX, maxIndexY, leftMaxIndexX,
                leftMaxIndexY, width, height);

  if (rightMax != 0)
    findHullTBB(S2, hull, maxIndexX, maxIndexY, rightX, rightY, rightMaxIndexX,
                rightMaxIndexY, width, height);
}

std::vector<int> quickHullTBB(std::vector<int> map, int width, int height,
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

  if (leftMax != 0)
    findHullTBB(S1, &hull, leftX, leftY, rightX, rightY, leftMaxIndexX,
                leftMaxIndexY, width, height);

  if (rightMax != 0)
    findHullTBB(S2, &hull, rightX, rightY, leftX, leftY, rightMaxIndexX,
                rightMaxIndexY, width, height);

  std::vector<int> result(hull);

  return result;
}

std::vector<int> convexQuickHullTBB(std::vector<int> map, int width,
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

  int hullsSize = hulls.size();

  tbb::task_scheduler_init init;

  tbb::parallel_for(tbb::blocked_range<int>(0, vert.size()),
                    [&](tbb::blocked_range<int> r) {
                      for (int i = r.begin(); i < r.end(); i++) {
                        hulls[i] = quickHullTBB(map, width, height, _marks[i],
                                                vert[i].first, vert[i].second);
                      }
                    });

  /*for (int i = 0; i < vertSize; i++) {
    hulls[i] = quickHullTBB(map, width, height, _marks[i], vert[i].first,
                            vert[i].second);
  }*/

  for (int i = 0; i < hullsSize; i++) {
    borders[i] = border++;
  }

  tbb::parallel_for(tbb::blocked_range<int>(0, hullsSize),
                    [&](tbb::blocked_range<int> r) {
                      for (int i = r.begin(); i < r.end(); i++) {
                        int hullSize = hulls[i].size();
                        for (int j = 0; j < hullSize; j++) {
                          result[hulls[i][j]] = borders[i];
                        }
                      }
                    });

  /*for (int i = 0; i < hullsSize; i++) {
    int hullSize = hulls[i].size();
    for (int j = 0; j < hullSize; j++) {
      result[hulls[i][j]] = borders[i];
    }
  }*/

  return result;
}
