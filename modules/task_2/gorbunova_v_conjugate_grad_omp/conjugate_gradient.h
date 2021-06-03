// Copyright 2021 Gorbunova Valeria
#ifndef MODULES_TASK_2_GORBUNOVA_V_CONJUGATE_GRAD_OMP_CONJUGATE_GRADIENT_H_
#define MODULES_TASK_2_GORBUNOVA_V_CONJUGATE_GRAD_OMP_CONJUGATE_GRADIENT_H_

double conj_grad(double** A, double* B, int S, int proc);
void MatrixMultiplicate(double** A, int S, double* Ax, double* x1);
double multiplicate(int S, double* x1, double* x2);
void Random_Matrix_A(int S, double** A);
void Random_Matrix_B(int S, double* B);
#endif  // MODULES_TASK_2_GORBUNOVA_V_CONJUGATE_GRAD_OMP_CONJUGATE_GRADIENT_H_
