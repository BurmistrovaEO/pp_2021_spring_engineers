// Copyright 2021 Kokh Vladislav
#include <omp.h>
#include <random>
#include <ctime>
#include "../../modules/task_2/kokh_v_hoar_sort_batcher_omp/hoar_sort_batcher.h"


int getDelta(int n, int p) {
    int dop = n / p;
    int a = dop + dop % 2;
    int dop_for_b = n - (p - 1) * a;
    int b = dop_for_b + dop_for_b % 2;
    int delta = (p - 1) * a + b - n;
    return delta;
}

double* randomFunc(int n, int p) {
    std::mt19937 gen(time(0));
    std::uniform_int_distribution <int> dist(-500, 500);
    if (p != -1) {
        n += getDelta(n, p);
    }
    double* mas = new double[n];
    for (int i = 0; i < n; ++i) {
        mas[i] = dist(gen);
    }
    return mas;
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

void splitter(double* mas, int start, int s1, int s2, int even) {
    int new_size2 = start + s1 + s2;
    int first = even;
    int second = start + s1 + even;
    int start_even = start + even;
    double* dop = new double[s1];
    for (int i = even; i < s1; i += 2) {
        dop[i] = mas[start + i];
    }
    while ((second < new_size2) && (first < s1)) {
        if (dop[first] >= mas[second]) {
            mas[start_even] = mas[second];
            second = second + 2;
        } else {
            mas[start_even] = dop[first];
            first = first + 2;
        }
        start_even = start_even + 2;
    }
    if (first < s1) {
        for (int i = first; i < s1; i += 2, start_even += 2) {
            mas[start_even] = dop[i];
        }
    } else {
        for (int i = second; i < new_size2; i += 2, start_even += 2) {
            mas[start_even] = mas[i];
        }
    }
}

void sort_parallel(double* mas, int n, int proc) {
    int  dop, a, dop_for_b, b, delta, new_n, m, proc_div2, last_value;
    omp_set_num_threads(proc);
    dop = n / proc;
    int int_max = 2147483647;
    a = dop % 2 + dop;
    dop_for_b = n - a * (proc - 1);
    b = dop_for_b + dop_for_b % 2;
    delta = (proc - 1) * a + b - n;
    new_n = n + delta;
    for (int i = n; i < new_n; ++i) {
        mas[i] = int_max;
    }
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < proc; ++i) {
            int right;
            int left = i * a;
            if (proc == i + 1) {
                right = delta + n - 1;
            } else {
                right = a - 1 + left;
            }
            sortFunc(left, right, mas);
        }
    }
    m = 1;
    proc_div2 = proc / 2;
    last_value = 0;
    while (proc_div2 > 0) {
#pragma omp parallel
        {
            int thread = omp_get_thread_num();
            if (thread < proc_div2 * 2) {
                if (thread % 2 == 0) {
                    int dop_param;
                    if (thread == proc - 2) {
                        dop_param = a * last_value + b;
                    } else {
                        dop_param = a * m;
                    }
                    splitter(mas, thread * a * m, a * m, dop_param, 0);
                } else {
                    int dop_param;
                    if (thread == proc - 1) {
                        dop_param = a * last_value + b;
                    } else {
                        dop_param = a * m;
                    }
                    splitter(mas, (thread - 1) * a * m, a * m, dop_param, 1);
                }
            }
        }
#pragma omp parallel
        {
            int thread = omp_get_thread_num();
            int start, f_value, dop_tmp;
            if (thread < 2 * proc_div2) {
                if (thread % 2 != 0) {
                    start = a * m * (thread - 1) + 3;
                    if (thread == proc - 1) {
                        f_value = new_n;
                    } else {
                        f_value = (thread + 1) * m * a;
                    }
                    for (int i = start; i < f_value - 1; i += 4) {
                        if (mas[i] > mas[i + 1]) {
                            dop_tmp = mas[i];
                            mas[i] = mas[i + 1];
                            mas[i + 1] = dop_tmp;
                        }
                    }
                } else {
                    start = thread * a * m + 1;
                    if (thread + 1 == proc - 1) {
                        f_value = new_n;
                    } else {
                        f_value = (thread + 2) * m * a;
                    }
                    for (int i = start; i < f_value - 1; i += 4) {
                        if (mas[i] > mas[i + 1]) {
                            dop_tmp = mas[i];
                            mas[i] = mas[i + 1];
                            mas[i + 1] = dop_tmp;
                        }
                    }
                }
            }
        }
        if (proc % 2 == 0) {
            last_value += m;
        }
        m = m * 2;
        proc = proc - proc_div2;
        proc_div2 = proc / 2;
    }
}
