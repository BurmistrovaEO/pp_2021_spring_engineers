// Copyright 2021 Kokh Vladislav
#ifndef MODULES_TASK_2_KOKH_V_HOAR_SORT_BATCHER_OMP_HOAR_SORT_BATCHER_H_
#define MODULES_TASK_2_KOKH_V_HOAR_SORT_BATCHER_OMP_HOAR_SORT_BATCHER_H_

double* randomFunc(int n, int p);
void sortFunc(int l, int r, double* mas);
bool check(double* mas, int n);
int getDelta(int n, int p);
void sort_parallel(double* mas, int n, int proc);
void splitter(double* mas, int start, int size1, int size2, int even);

#endif     // MODULES_TASK_2_KOKH_V_HOAR_SORT_BATCHER_OMP_HOAR_SORT_BATCHER_H_
