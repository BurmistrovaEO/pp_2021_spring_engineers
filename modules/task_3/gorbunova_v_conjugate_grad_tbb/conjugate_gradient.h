// Copyright 2021 Gorbunova Valeria
#ifndef MODULES_TASK_3_GORBUNOVA_V_CONJUGATE_GRAD_TBB_CONJUGATE_GRADIENT_H_
#define MODULES_TASK_3_GORBUNOVA_V_CONJUGATE_GRAD_TBB_CONJUGATE_GRADIENT_H_

double conj_gradP(double** A, double* B, int S, int proc);
void MatrixMultiplicateP(double** A, int S, double* Ax, double* x1, int grainsize);
double multiplicate(int S, double* x1, double* x2);
void Random_Matrix_A(int S, double** A);
void Random_Matrix_B(int S, double* B);
#endif  // MODULES_TASK_3_GORBUNOVA_V_CONJUGATE_GRAD_TBB_CONJUGATE_GRADIENT_H_
