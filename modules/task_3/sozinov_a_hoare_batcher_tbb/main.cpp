// Copyright 2021 Sozinov Alex

#include <gtest/gtest.h>
#include <tbb/tbb.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "./hoare_batcher_tbb.h"

TEST(GetRandomVector, get_vector) {
    ASSERT_NO_THROW(GetRandomVector(110));
}

TEST(GetRandomVector, vector_have_correct_size) {
    std::vector<double> vect = GetRandomVector(100);
    ASSERT_EQ(static_cast<int>(vect.size()), 100);
}

TEST(Sort, can_sort_one_thread) {
    std::vector<double> vect = GetRandomVector(100);
    ASSERT_NO_THROW(ParSort(&vect, 1));
}

TEST(Sort, can_sort_four_thread) {
    std::vector<double> vect = GetRandomVector(100);
    ASSERT_NO_THROW(ParSort(&vect, 4));
}

TEST(Sort, correct_sort_10000) {
    std::vector<double> vect = GetRandomVector(100);
    std::vector<double> copy = vect;
    std::vector<double> par_vect = vect;
    std::sort(copy.begin(), copy.end());
    tbb::tick_count seq_start = tbb::tick_count::now();
    SeqSort(&vect);
    tbb::tick_count seq_end = tbb::tick_count::now();

    tbb::tick_count par_start = tbb::tick_count::now();
    ParSort(&par_vect, 1);
    tbb::tick_count par_end = tbb::tick_count::now();

    std::cout << "vector size = " << copy.size() << std::endl;
    std::cout << "Seq: time: " << (seq_end - seq_start).seconds()
              << "\tresult vector size: " << vect.size() << std::endl;
    std::cout << "Par: time: " << (par_end - par_start).seconds()
              << "\tresult vector size: " << par_vect.size() << std::endl;

    ASSERT_EQ(copy, vect);
    ASSERT_EQ(copy, par_vect);
}

TEST(Sort, correct_sort_900000) {
    std::vector<double> vect = GetRandomVector(900);
    std::vector<double> copy = vect;
    std::vector<double> par_vect = vect;
    std::sort(copy.begin(), copy.end());
    tbb::tick_count seq_start = tbb::tick_count::now();
    SeqSort(&vect);
    tbb::tick_count seq_end = tbb::tick_count::now();

    tbb::tick_count par_start = tbb::tick_count::now();
    ParSort(&par_vect, 4);
    tbb::tick_count par_end = tbb::tick_count::now();

    std::cout << "vector size = " << copy.size() << std::endl;
    std::cout << "Seq: time: " << (seq_end - seq_start).seconds()
              << "\tresult vector size: " << vect.size() << std::endl;
    std::cout << "Par: time: " << (par_end - par_start).seconds()
              << "\tresult vector size: " << par_vect.size() << std::endl;

    ASSERT_EQ(copy, vect);
    ASSERT_EQ(copy, par_vect);
}

TEST(Sort, correct_sort_1000000) {
    std::vector<double> vect = GetRandomVector(100);
    std::vector<double> copy = vect;
    std::vector<double> par_vect = vect;
    std::sort(copy.begin(), copy.end());
    tbb::tick_count seq_start = tbb::tick_count::now();
    SeqSort(&vect);
    tbb::tick_count seq_end = tbb::tick_count::now();

    tbb::tick_count par_start = tbb::tick_count::now();
    ParSort(&par_vect, 4);
    tbb::tick_count par_end = tbb::tick_count::now();

    std::cout << "vector size = " << copy.size() << std::endl;
    std::cout << "Seq: time: " << (seq_end - seq_start).seconds()
              << "\tresult vector size: " << vect.size() << std::endl;
    std::cout << "Par: time: " << (par_end - par_start).seconds()
              << "\tresult vector size: " << par_vect.size() << std::endl;

    ASSERT_EQ(copy, vect);
    ASSERT_EQ(copy, par_vect);
}

TEST(Sort, correct_sort_3000000) {
    std::vector<double> vect = GetRandomVector(300);
    std::vector<double> copy = vect;
    std::vector<double> par_vect = vect;
    std::sort(copy.begin(), copy.end());
    tbb::tick_count seq_start = tbb::tick_count::now();
    SeqSort(&vect);
    tbb::tick_count seq_end = tbb::tick_count::now();

    tbb::tick_count par_start = tbb::tick_count::now();
    ParSort(&par_vect, 2);
    tbb::tick_count par_end = tbb::tick_count::now();

    std::cout << "vector size = " << copy.size() << std::endl;
    std::cout << "Seq: time: " << (seq_end - seq_start).seconds()
              << "\tresult vector size: " << vect.size() << std::endl;
    std::cout << "Par: time: " << (par_end - par_start).seconds()
              << "\tresult vector size: " << par_vect.size() << std::endl;

    ASSERT_EQ(copy, vect);
    ASSERT_EQ(copy, par_vect);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
