// Copyright 2021 Ekaterina Burmistrova
#include <tbb/tbb.h>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <functional>
#include "./../../modules/task_3/bourmistrova_e_trapez_integr_tbb/tr_int_tbb.h"


double CheckCoeff(double i, int s_pr) {
    double qju = 1.00;
    if (i = 0 || i == s_pr)
        qju *= 0.5;
    return qju;
}

double OneDimIntegr(const std::vector<std::pair<int, int>>& bord,
    std::function<double(double, double, double)> f,
    int set_prec, int pr1, double pr2, double s, int b) {
    //// tbb::task_scheduler_init(6);  // automatic
    double tr_sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, set_prec),
        0.0,
        [&](tbb::blocked_range<int> r, double running_total){
            for (pr1 = r.begin(); pr1 < r.end(); ++pr1) {
                pr2 = bord[b].first + pr1 * s;
                running_total += ((f(pr2, 1, 1) + f(pr2 + s, 1, 1)) / 2) * s;
            }
            return running_total;
        }, std::plus<double>());
    return tr_sum;
}
double SolveParallel(const std::vector<std::pair<int, int>>& bord,
    std::function<double(double, double, double)> f) {
    //// tbb::task_scheduler_init(6);  // automatic
    double tr_sum = 0;
    double m = 1;
    int i = 0, j = 0, k = 0;
    double q0 = 0, q1 = 0, q2 = 0;
    double x = 0, x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0;
    int set_precision = 1000;
    double s = (bord[0].second -
        bord[0].first) / static_cast<double>(set_precision);
    if (bord.size() == 1) {  // ONE DIMENSION
        // std::cout << "In OneDim ";
        tr_sum = OneDimIntegr(bord, f, set_precision, i, x, s, 0);
        // std::cout << tr_sum;
    } else if (bord.size() == 2) {  // TWO DIMENSIONS
        double s2 = (bord[1].second -
            bord[1].first) / static_cast<double>(set_precision);
                for (i = 0; i < set_precision; ++i) {
                    tr_sum += tbb::parallel_reduce(
                        tbb::blocked_range<int>(0, set_precision), 0.0,
                        [&](tbb::blocked_range<int> r, double running_total) {
                            for (j = r.begin(); j < r.end(); ++j) {
                                x = bord[0].first + i * s;
                                x1 = x + i * s;
                                x2 = bord[1].first + j * s2;
                                x3 = x2 + j * s2;
                                // if (i == 0 || i == set_precision)
                                //    m *= 0.5;
                                // if (j == 0 || j == set_precision)
                                //   m *= 0.5;
                                running_total += (f(x, x2, 1) + f(x, x3, 1) +
                                    f(x1, x2, 1) + f(x1, x3, 1));
                                // m = 1;
                            }
                            return running_total;
                        }, std::plus<double>());
            }
        tr_sum = ((s * s2)/4) * tr_sum;
    } else {  // THREE DIMENSIONS
        double s2 = (bord[1].second -
            bord[1].first) / static_cast<double>(set_precision);
        double s3 = (bord[2].second -
            bord[2].first) / static_cast<double>(set_precision);
        double q0 = 1, q1 = 1, q2 = 1;
        double tmp = 0;
        tr_sum += tbb::parallel_reduce(tbb::blocked_range2d<int>(0,
            set_precision, 0, set_precision), 0.0,
            [&](tbb::blocked_range2d<int, int> r, double running_total) {
                for (i = r.rows().begin(); i < r.rows().end(); ++i) {
        x = bord[1].first + i * s2;
     x1 = x + i * s2;
      for (j = r.cols().begin(); j < r.cols().end(); ++j) {
           x2 = bord[2].first + j * s3;
      x3 = x2 + j * s3;
             running_total += (f(x2, x, 1) + f(x2, x1, 1) +
                 f(x3, x, 1) + f(x3, x1, 1))/2;
                  }
          } return running_total;
         }, std::plus<double>());
        tr_sum += tbb::parallel_reduce(tbb::blocked_range<int>(0,
        set_precision), 0.0,
        [&](tbb::blocked_range<int> r, double running_total) {
        for (i = r.begin(); i < r.end(); ++i) {
               x = bord[0].first + i * s;
                 x1 = x2 + i * s;
               running_total += (f(1, 1, x) + f(1, 1, x1))/2;
           } return running_total;
           }, std::plus<double>());
           tr_sum = ((s * s2)/4) * tr_sum;
    }
    return tr_sum;
}
