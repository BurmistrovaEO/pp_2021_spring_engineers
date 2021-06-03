// Copyright 2021 Pestrikov Ilya
#include <omp.h>
#include <cstdlib>
#include <vector>
#include <random>
#include <iostream>
#include <cmath>
#include <ctime>
#include "./pestrikov_i_graham_omp.h"

std::vector<Point>* graham(std::vector<Point>* points, int startIndex,
    int endIndex) {
    double start = omp_get_wtime();

    Point leftBottomPoint = findLeftBottomPoint(*points, startIndex, endIndex);
    moveCoordByXY(points, -leftBottomPoint.x, -leftBottomPoint.y, startIndex,
        endIndex);

    sort((*points).begin() + startIndex, (*points).begin() + endIndex,
        [&](Point left, Point right) {
        return multiplyVectors(left, right) > 0 ||
            (multiplyVectors(left, right) == 0 &&
            left.x * left.x + left.y * left.y <
            right.x * right.x + right.y * right.y);
        });

    std::vector<Point>* hull = new std::vector<Point>();
    for (int i = startIndex; i < endIndex; ++i) {
        while (hull->size() >= 2 &&
            multiplyVectors(
                substractPoints((*points)[i], (*hull)[hull->size() - 1]),
                substractPoints(
                    (*hull)[hull->size() - 2],
                    (*hull)[hull->size() - 1])) <= 0) {
            hull->pop_back();
        }
        hull->push_back((*points)[i]);
    }

    moveCoordByXY(hull, leftBottomPoint.x, leftBottomPoint.y, 0,
        hull->size());

    double end = omp_get_wtime();
    if (omp_in_parallel() == 0) {
        std::cout << "Graham sequential time: " << (end - start) * 1000 <<
            " ms" << std::endl;
    }

    return hull;
}

std::vector<Point>* grahamOmp(std::vector<Point>* points,
    int threadsCount) {
    double start = omp_get_wtime();

    int step = (*points).size() / threadsCount;

    std::vector<Point>* combinedHull = new std::vector<Point>();
    #pragma omp parallel num_threads(threadsCount)
    {
        std::vector<Point>* threadHull;
        int currentThreadNumber = omp_get_thread_num();

        int startIndex = currentThreadNumber * step;
        int endIndex = (currentThreadNumber + 1) * step;
        if (currentThreadNumber == threadsCount - 1) {
            endIndex = (*points).size();
        }
        threadHull = graham(points, startIndex, endIndex);

        #pragma omp critical
        {
            for (auto point : *threadHull) {
                combinedHull->push_back(point);
            }
        }
    }

    std::vector<Point>* result = graham(combinedHull, 0,
        combinedHull->size());

    double end = omp_get_wtime();
    std::cout << "Graham omp time (threadsCount: " << threadsCount << "): " <<
        (end - start) * 1000 << " ms" << std::endl;

    return result;
}

Point findLeftBottomPoint(const std::vector<Point>& points, int startIndex,
    int endIndex) {
    Point leftBottomPoint = points[startIndex];

    for (int i = startIndex; i < endIndex; ++i) {
        if (points[i].x < leftBottomPoint.x ||
            (points[i].x == leftBottomPoint.x &&
            points[i].y < leftBottomPoint.y)) {
            leftBottomPoint = points[i];
        }
    }

    return leftBottomPoint;
}

void moveCoordByXY(std::vector<Point>* points, int x, int y,
    int startIndex, int endIndex) {
    for (int i = startIndex; i < endIndex; ++i) {
        (*points)[i].x += x;
        (*points)[i].y += y;
    }
}

int multiplyVectors(const Point& point1, const Point& point2) {
    return point1.x * point2.y - point1.y * point2.x;
}

Point substractPoints(const Point& point1, const Point& point2) {
    return { point1.x - point2.x, point1.y - point2.y };
}

bool isPointInsidePolygon(const std::vector<Point>& polygon,
    const Point& point) {
    bool result = false;
    for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        if (((polygon[i].y >= point.y) != (polygon[j].y >= point.y)) &&
            (point.x <=
                static_cast<double>(polygon[j].x - polygon[i].x) *
                (point.y - polygon[i].y) / (polygon[j].y - polygon[i].y) +
                polygon[i].x)) {
            result = !result;
        }
    }
    return result;
}

bool isPointBetweenPoints(Point a, Point b, Point c) {
    int crossProduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y);

    if (std::abs(crossProduct) != 0) {
        return false;
    }

    int dotProduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y) * (b.y - a.y);
    if (dotProduct < 0) {
        return false;
    }

    int squaredLength = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
    if (dotProduct > squaredLength) {
        return false;
    }

    return true;
}

bool isPointOnPolygonBorder(const std::vector<Point>& polygon,
    const Point& point) {
    for (auto polygonPoint : polygon) {
        if (point.x == polygonPoint.x && point.y == polygonPoint.y) {
            return true;
        }
    }

    for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
        if (isPointBetweenPoints(polygon[i], polygon[j], point)) {
            return true;
        }
    }

    return false;
}

bool isGrahamCorrect(const std::vector<Point>& polygon,
    const std::vector<Point>& points) {
    for (auto point : points) {
        if (!isPointInsidePolygon(polygon, point) &&
            !isPointOnPolygonBorder(polygon, point)) {
            return false;
        }
    }

    return true;
}

int getRandomInt(int min, int max) {
    return (std::rand() % (max - min + 1)) + min;
}



std::vector<Point>* generateRandomPointsVector(int size, int min, int max) {
    std::srand(time(0));
    std::vector<Point>* result = new std::vector<Point>();

    for (int i = 0; i < size; ++i) {
        result->push_back({ getRandomInt(min, max), getRandomInt(min, max) });
    }

    return result;
}
