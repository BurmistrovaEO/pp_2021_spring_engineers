// Copyright 2021 Pestrikov Ilya
#pragma once

#include <vector>

struct Point {
    int x;
    int y;
};

std::vector<Point>* graham(std::vector<Point>* points, int startIndex,
    int endIndex);

std::vector<Point>* grahamTbb(std::vector<Point>* points);

Point findLeftBottomPoint(const std::vector<Point>& points, int startIndex,
    int endIndex);

void moveCoordByXY(std::vector<Point>* points, int x, int y,
    int startIndex, int endIndex);

int multiplyVectors(const Point& point1, const Point& point2);

Point substractPoints(const Point& point1, const Point& point2);

bool isPointInsidePolygon(const std::vector<Point>& polygon,
    const Point& point);

bool isPointOnPolygonBorder(const std::vector<Point>& polygon,
    const Point& point);

bool isGrahamCorrect(const std::vector<Point>& polygon,
    const std::vector<Point>& points);

int getRandomInt(int min, int max);

std::vector<Point>* generateRandomPointsVector(int size, int min, int max);
