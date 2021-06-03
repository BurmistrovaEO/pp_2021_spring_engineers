// Copyright 2021 Ekaterina Burmistrova
#include <omp.h>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <functional>
#include "./../../modules/task_2/bourmistrova_e_trapez_integr_omp/trapezoid_integral_omp.h"


double CheckCoeff(double i, int s_pr) {
    double qju = 1.00;
    if (i = 0 || i == s_pr)
        qju *= 0.5;
    return qju;
}

double OneDimIntegr(const std::vector<std::pair<int, int>>& bord,
    std::function<double(double, double, double)> f,
    int set_prec, int pr1, double pr2, double s, int b) {
    double tr_sum = 0;
    omp_set_num_threads(6);
#pragma omp parallel for private(pr1, pr2) reduction(+:tr_sum)
    for (pr1 = 0; pr1 < set_prec; pr1++) {
        pr2 = bord[b].first + pr1 * s;
        tr_sum += ((f(pr2, 1, 1) + f(pr2 + s, 1, 1)) / 2) * s;
    }
    return tr_sum;
}
double SolveParallel(const std::vector<std::pair<int, int>>& bord,
    std::function<double(double, double, double)> f) {
    omp_set_num_threads(4);
    double tr_sum = 0;
    double m = 1;
    int i = 0, j = 0, k = 0;
    double q0 = 0, q1 = 0, q2 = 0;
    double x = 0, x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0;
    int set_precision = 1000;
    double s = (bord[0].second -
        bord[0].first) / static_cast<double>(set_precision);
    if (bord.size() == 1) {  // ONE DIMENSION
        tr_sum = OneDimIntegr(bord, f, set_precision, i, x, s, 0);
    } else if (bord.size() == 2) {  // TWO DIMENSIONS
        double s2 = (bord[1].second -
            bord[1].first) / static_cast<double>(set_precision);
#pragma omp parallel for reduction(+:tr_sum)
        for (i = 0; i < set_precision + 1; i++)
            for (j = 0; j < set_precision + 1; j++) {
                x = bord[0].first + i * s;
                x1 = x + i * s;
                x2 = bord[1].first + j * s2;
                x3 = x2 + j * s2;
                tr_sum += (f(x, x2, 1) + f(x, x3, 1) +
                    f(x1, x2, 1) + f(x1, x3, 1));
            }
        tr_sum = ((s * s2) / 4) * tr_sum;
    } else {  // THREE DIMENSIONS
        double s2 = (bord[1].second -
            bord[1].first) / static_cast<double>(set_precision);
        double s3 = (bord[2].second -
            bord[2].first) / static_cast<double>(set_precision);
        double q0 = 1, q1 = 1, q2 = 1;
        double tmp = 0;
#pragma omp parallel for reduction(+:tr_sum)
        for (i = 0; i < set_precision; i++) {
            x = bord[2].first + i * s3;
            tr_sum += ((f(x, 1, 1) + f(x + s3, 1, 1)) / 2)*s3;
        }
#pragma omp parallel for reduction(+:tmp)
        for (i = 0; i < set_precision; i++) {
            x = bord[1].first + i * s2;
            tmp += ((f(1, x, 1) + f(1, x + s2, 1)) / 2)*s2 * tr_sum;
        }
#pragma omp parallel for reduction(+:tr_sum)
        for (i = 0; i < set_precision; i++) {
            x = bord[0].first + i * s;
            tr_sum += ((f(1, 1, x) + f(1, 1, x + s3)) / 2)*s*tmp;
        }
    }
    return tr_sum;
}
double SolveParallelSum(const std::vector<std::pair<int, int>>& bord,
    std::function<double(double, double, double)> f) {
    omp_set_num_threads(4);
    double tr_sum = 0;
    double m = 1;
    int i = 0, j = 0, k = 0;
    double q0 = 0, q1 = 0, q2 = 0;
    double x = 0, x1 = 0, x2 = 0, x3 = 0;
    int set_precision = 1000;
    double s = (bord[0].second -
        bord[0].first) / static_cast<double>(set_precision);
    if (bord.size() == 2) {  // TWO DIMENSIONS
        double s2 = (bord[1].second -
            bord[1].first) / static_cast<double>(set_precision);
        for (i = 0; i < set_precision + 1; i++)
#pragma omp parallel for schedule(static) reduction(+:tr_sum)
            for (j = 0; j < set_precision + 1; j++) {
                x = bord[0].first + i * s;
                x1 = x + i * s;
                x2 = bord[1].first + j * s2;
                x3 = x2 + j * s2;
                if (i == 0 || i == set_precision)
                    m *= 0.5;
                if (j == 0 || j == set_precision)
                    m *= 0.5;
                tr_sum += m * (f(x, x2, 1) + f(x1, x3, 1) +
                    f(x1, x2, 1) + f(x1, x3, 1));
                m = 1;
            }
        tr_sum = ((s * s2) / 7.5) * tr_sum;
    } else {  // THREE DIMENSIONS
        double s2 = (bord[1].second -
            bord[1].first) / static_cast<double>(set_precision);
        double s3 = (bord[2].second -
            bord[2].first) / static_cast<double>(set_precision);
        double m = 1;
        double q0 = 1, q1 = 1, q2 = 1;
        double tmp = 0;
#pragma omp parallel for reduction(+:tr_sum)
        for (i = 0; i < set_precision; i++) {
            x = bord[2].first + i * s3;
            tr_sum += ((f(x, 1, 1) + f(x + s3, 1, 1)) / 2) * s3;
        }
#pragma omp parallel for reduction(+:tr_sum)
        for (i = 0; i < set_precision; i++) {
            x = bord[1].first + i * s2;
            tr_sum += ((f(1, x, 1) + f(1, x + s2, 1)) / 2) * s2;
        }
#pragma omp parallel for reduction(+:tr_sum)
        for (i = 0; i < set_precision; i++) {
            x = bord[0].first + i * s;
            tr_sum += ((f(1, 1, x) + f(1, 1, x + s3)) / 2) * s;
        }
    }
    return tr_sum;
}


