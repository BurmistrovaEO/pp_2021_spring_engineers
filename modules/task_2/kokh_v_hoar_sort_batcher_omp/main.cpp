// Copyright 2021 Kokh Vladislav
#include <gtest/gtest.h>
#include "./hoar_sort_batcher.h"

TEST(hoar_sort_batcher, simple_check_5000) {
    int size = 5000;
    double* mas = randomFunc(size, 10);
    sort_parallel(mas, size, 10);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}

TEST(hoar_sort_batcher, check_otr) {
    int size = 10000;
    double* mas = randomFunc(size, 8);
    mas[7] = -100;
    sort_parallel(mas, size, 8);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}
TEST(hoar_sort_batcher, simple_check_10000) {
    int size = 10000;
    double* mas = randomFunc(size, 5);
    sort_parallel(mas, size, 5);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}
TEST(hoar_sort_batcher, check_all_otr) {
    int size = 5000;
    double* mas = randomFunc(size, 1);
    sort_parallel(mas, size, 1);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}
TEST(hoar_sort_batcher, simple_check_all_zero) {
    int size = 5000;
    double* mas = randomFunc(size, 5);
    sort_parallel(mas, size, 5);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
    double* mas2 = randomFunc(size, -1);
    sortFunc(0, size - 1, mas2);
    ASSERT_EQ(true, check(mas2, size));
    delete[] mas2;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
