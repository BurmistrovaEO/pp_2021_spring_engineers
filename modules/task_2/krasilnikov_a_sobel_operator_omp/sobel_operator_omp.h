// Copyright 2021 Krasilnikov Alexey
#ifndef MODULES_TASK_2_KRASILNIKOV_A_SOBEL_OPERATOR_OMP_SOBEL_OPERATOR_OMP_H_
#define MODULES_TASK_2_KRASILNIKOV_A_SOBEL_OPERATOR_OMP_SOBEL_OPERATOR_OMP_H_


void SobelOperator_parallel(unsigned char* image, int height, int weight);

void SobelOperator(unsigned char* image, int height, int weight);

#endif  // MODULES_TASK_2_KRASILNIKOV_A_SOBEL_OPERATOR_OMP_SOBEL_OPERATOR_OMP_H_
