// Copyright 2021 Kokh Vladislav


#include <tbb/blocked_range.h>
#include <tbb/task.h>
#include <tbb/tbb.h>
#include <random>
#include <algorithm>
#include <ctime>
#include <cmath>

#include "../../modules/task_3/kokh_v_hoar_sort_batcher_tbb/hoar_sort_batcher.h"


double* randomFunc(int n) {
    std::mt19937 gen(time(0));
    std::uniform_int_distribution <int> dist(-500, 500);
    double* mas = new double[n];
    for (int i = 0; i < n; ++i) {
        mas[i] = dist(gen);
    }
    return mas;
}

bool check(double* mas, int n) {
    bool flag = true;
    for (int i = 0; i < n - 1; i++) {
        if (mas[i] > mas[i + 1]) {
            flag = false;
            break;
        }
    }
    return flag;
}

int get_part(int n, int count_of_threads) {
    int part = ceil(n / count_of_threads);
    return part;
}

void tbb_sort_main(double* mas, int n, int count_of_threads) {
    int part = get_part(n, count_of_threads);
    double* dop = new double[n];
    tbb_sorter& sr = *new(tbb::task::allocate_root())tbb_sorter(mas,
    dop, n, part);
    tbb::task::spawn_root_and_wait(sr);
    delete[] dop;
}


void sortFunc(int l, int r, double* mas) {
    int l_frame, r_frame;
    double razr_el, dop;
    while (r > l) {
        l_frame = l;
        r_frame = r;
        razr_el = mas[(l + r) / 2];
        while (l_frame <= r_frame) {
            while (mas[l_frame] < razr_el) {
                l_frame++;
            }
            while (mas[r_frame] > razr_el) {
                r_frame--;
            }
            if (l_frame <= r_frame) {
                dop = mas[l_frame];
                mas[l_frame] = mas[r_frame];
                mas[r_frame] = dop;
                l_frame++;
                r_frame--;
            }
        }
        if (2 * l_frame > l + r) {
            sortFunc(l_frame, r, mas);
            r = l_frame - 1;
        } else {
            sortFunc(l, l_frame - 1, mas);
            l = l_frame;
        }
    }
}

tbb::task* splitter::execute() {
    int first = even;
    int second = even;
    int start_even = even;
    double* mas_second = mas + size1;
    for (int i = even; i < size1; i += 2) {
        dop[i] = mas[i];
    }
    while ((second < size2) && (first < size1)) {
        if (dop[first] >= mas_second[second]) {
            mas[start_even] = mas_second[second];
            second = second + 2;
        } else {
            mas[start_even] = dop[first];
            first = first + 2;
        }
        start_even = start_even + 2;
    }
    if (first != size1) {
        for (int i = first; i < size1; i += 2, start_even += 2) {
            mas[start_even] = dop[i];
        }
    } else {
        for (int i = second; i < size2; i += 2, start_even += 2) {
            mas[start_even] = mas_second[i];
        }
    }
    return NULL;
}


void compare::operator()(const tbb::blocked_range<int>& range_value) const {
    int begin, end, double_i, double_i_min;
    begin = range_value.begin();
    end = range_value.end();
    for (int i = begin; i < end; i++) {
        double_i = 2 * i;
        double_i_min = 2 * i - 1;
        if (mas[double_i] < mas[double_i_min]) {
            std::swap(mas[double_i_min], mas[double_i]);
        }
    }
}


tbb::task* tbb_sorter::execute() {
    if (n > part) {
        int hauf_size, s, de;
        hauf_size = n / 2;
        s = hauf_size + hauf_size % 2;
        tbb_sorter& sr1 = *new(allocate_child())tbb_sorter(mas, dop, s, part);
        tbb_sorter& sr2 = *new(allocate_child())tbb_sorter(mas + s,
        dop + s, n - s, part);
        set_ref_count(3);
        spawn(sr1);
        spawn_and_wait_for_all(sr2);
        splitter& sp1 = *new(allocate_child()) splitter(mas, dop, s, n - s, 0);
        splitter& sp2 = *new(allocate_child()) splitter(mas, dop, s, n - s, 1);
        set_ref_count(3);
        spawn(sp1);
        spawn_and_wait_for_all(sp2);
        de = (n + 1) / 2;
        tbb::parallel_for(tbb::blocked_range<int>(1, de), compare(mas));
    } else { sortFunc(0, n - 1, mas); }
    return NULL;
}
