// Copyright 2021 Solomakhin Sergey
#include <omp.h>
#include "gtest/gtest.h"



#include "./hoar_sort_simple_merge.h"

TEST(randomgen, run) {
    ASSERT_NO_THROW(random_gen(100));
}

TEST(randomgen, wrongsize) {
    ASSERT_ANY_THROW(random_gen(-100));
}

TEST(hoarsort, runwithoutexceptionswithsize10) {
    auto vect = random_gen(10);
    ASSERT_NO_THROW(hoar_sort_omp(0, 9, &vect));
}

TEST(hoarsort, runwithoutexceptionswithsize_100) {
    auto vect = random_gen(100);
    ASSERT_NO_THROW(hoar_sort_omp(0, 99, &vect));
}

TEST(hoarsort, runwithoutexceptionswithsize1000) {
    auto vect = random_gen(1000);
    ASSERT_NO_THROW(hoar_sort_omp(0, 999, &vect));
}

TEST(hoarsort, runwithsize10) {
    auto vect = random_gen(10);
    auto vect2 = random_gen(10);
    double start = omp_get_wtime();
    hoar_sort(0, 9, &vect);
    double end = omp_get_wtime();
    std::cout << "Seq sort time -> " << end - start << std::endl;
    start = omp_get_wtime();
    hoar_sort_omp(0, 9, &vect2);
    end = omp_get_wtime();
    std::cout << "OMP sort time -> " << end - start << std::endl;
    ASSERT_EQ(1, checker(&vect2));
}

TEST(hoarsort, runwithsize100) {
    auto vect = random_gen(100);
    auto vect2 = random_gen(100);
    double start = omp_get_wtime();
    hoar_sort(0, 99, &vect);
    double end = omp_get_wtime();
    std::cout << "Seq sort time -> " << end - start << std::endl;
    start = omp_get_wtime();
    hoar_sort_omp(0, 99, &vect2);
    end = omp_get_wtime();
    std::cout << "OMP sort time -> " << end - start << std::endl;
    ASSERT_EQ(1, checker(&vect2));
}

TEST(hoarsort, runrandomsize1000) {
    auto vect = random_gen(1000);
    auto vect2 = random_gen(1000);
    double start = omp_get_wtime();
    hoar_sort(0, 999, &vect);
    double end = omp_get_wtime();
    std::cout << "Seq sort time -> " << end - start << std::endl;
    start = omp_get_wtime();
    hoar_sort_omp(0, 999, &vect2);
    end = omp_get_wtime();
    std::cout << "OMP sort time -> " << end - start << std::endl;
    ASSERT_EQ(1, checker(&vect2));
}

TEST(hoarsort, runrandomsize10000) {
    auto vect = random_gen(10000);
    auto vect2 = random_gen(10000);
    double start = omp_get_wtime();
    hoar_sort(0, 9999, &vect);
    double end = omp_get_wtime();
    std::cout << "Seq sort time -> " << end - start << std::endl;
    start = omp_get_wtime();
    hoar_sort_omp(0, 9999, &vect2);
    end = omp_get_wtime();
    std::cout << "OMP sort time -> " << end - start << std::endl;
    ASSERT_EQ(1, checker(&vect2));
}

TEST(hoarsort, runrandomsize100000) {
    auto vect = random_gen(100000);
    auto vect2 = random_gen(100000);
    double start = omp_get_wtime();
    hoar_sort(0, 99999, &vect);
    double end = omp_get_wtime();
    std::cout << "Seq sort time -> " << end - start << std::endl;
    start = omp_get_wtime();
    hoar_sort_omp(0, 99999, &vect2);
    end = omp_get_wtime();
    std::cout << "OMP sort time -> " << end - start << std::endl;
    ASSERT_EQ(1, checker(&vect2));
}

TEST(merging, runwithsize10tosize20) {
    auto vect_a = random_gen(10);
    auto vect_b = random_gen(10);
    hoar_sort_omp(0, 9, &vect_a);
    hoar_sort_omp(0, 9, &vect_b);
    auto merged = merge(vect_a, vect_b);
    for (int i = 0; i < static_cast<int>(vect_a.size()) - 1; ++i) {
        ASSERT_TRUE(merged[i] <= merged[i + 1]);
    }
}

TEST(merging, runwithsize100tosize200) {
    auto vect_a = random_gen(100);
    auto vect_b = random_gen(100);
    hoar_sort_omp(0, 99, &vect_a);
    hoar_sort_omp(0, 99, &vect_b);
    auto merged = merge(vect_a, vect_b);
    for (int i = 0; i < static_cast<int>(vect_a.size()) - 1; ++i) {
        ASSERT_TRUE(merged[i] <= merged[i + 1]);
    }
}

TEST(merging, runwithsize1000tosize2000) {
    auto vect_a = random_gen(1000);
    auto vect_b = random_gen(1000);
    hoar_sort_omp(0, 999, &vect_a);
    hoar_sort_omp(0, 999, &vect_b);
    auto merged = merge(vect_a, vect_b);
    for (int i = 0; i < static_cast<int>(vect_a.size()) - 1; ++i) {
        ASSERT_TRUE(merged[i] <= merged[i + 1]);
    }
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
