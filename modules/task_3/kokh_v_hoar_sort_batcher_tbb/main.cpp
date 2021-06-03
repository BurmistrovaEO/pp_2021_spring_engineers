// Copyright 2021 Kokh Vladislav
#include <gtest/gtest.h>
#include "./hoar_sort_batcher.h"


TEST(hoar_sort_batcher, simple_check_10000_44) {
    int size = 100000;
    double* mas = randomFunc(size);
    tbb_sort_main(mas, size, 4);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}

TEST(hoar_sort_batcher, simple_check_5000_5) {
    int size = 5000;
    double* mas = randomFunc(size);
    tbb_sort_main(mas, size, 5);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}

TEST(hoar_sort_batcher, check_otr) {
    int size = 10000;
    double* mas = randomFunc(size);
    mas[7] = -100;
    tbb_sort_main(mas, size, 8);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}

TEST(hoar_sort_batcher, simple_check_10000) {
    int size = 10000;
    double* mas = randomFunc(size);
    tbb_sort_main(mas, size, 6);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}

TEST(hoar_sort_batcher, check_all_otr) {
    int size = 5000;
    double* mas = randomFunc(size);
    tbb_sort_main(mas, size, 1);
    ASSERT_EQ(true, check(mas, size));
    delete[] mas;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
