// Copyright 2021 Krasilnikov Alexey
#ifndef MODULES_TASK_3_KRASILNIKOV_A_SOBEL_OPERATOR_TBB_SOBEL_OPERATOR_TBB_H_
#define MODULES_TASK_3_KRASILNIKOV_A_SOBEL_OPERATOR_TBB_SOBEL_OPERATOR_TBB_H_


void SobelOperator_parallel(unsigned char* image, int height, int weight);

void SobelOperator(unsigned char* image, int height, int weight);

#endif  // MODULES_TASK_3_KRASILNIKOV_A_SOBEL_OPERATOR_TBB_SOBEL_OPERATOR_TBB_H_
