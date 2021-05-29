// Copyright 2021 Kokh Vladislav
#ifndef MODULES_TASK_3_KOKH_V_HOAR_SORT_BATCHER_TBB_HOAR_SORT_BATCHER_H_
#define MODULES_TASK_3_KOKH_V_HOAR_SORT_BATCHER_TBB_HOAR_SORT_BATCHER_H_
#include <tbb/blocked_range.h>
#include <tbb/task.h>


int get_part(int n, int count_of_threads);

class compare {
 private:
    double* mas;
 public:
    explicit compare(double* m) {
        mas = m;
    }    void operator()(const tbb::blocked_range<int>& range_value) const;
};

class splitter :public tbb::task {
 private:
    double* mas;
    double* dop;
    int size1;
    int size2;
    int even;
 public:
    splitter(double* m, double* d, int s1, int s2, int e) {
        mas = m;
        dop = d;
        size1 = s1;
        size2 = s2;
        even = e;
    }    tbb::task* execute();
};



void sortFunc(int l, int r, double* mas);
double* randomFunc(int n);
bool check(double* mas, int n);

class tbb_sorter :public tbb::task {
 private:
    double* mas;
    double* dop;
    int n;
    int part;
 public:
    tbb_sorter(double* m, double* d, int s, int p) {
        mas = m;
        dop = d;
        n = s;
        part = p;
    }
    tbb::task* execute();
};

void tbb_sort_main(double* mas, int n, int count_of_threads);


#endif  // MODULES_TASK_3_KOKH_V_HOAR_SORT_BATCHER_TBB_HOAR_SORT_BATCHER_H_
