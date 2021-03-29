#include <arm_neon.h>

#ifndef UTILS_H
#define UTILS_H

double average(uint64_t nums[], uint64_t len);

uint64_t compare_doubles(double a, double b, double errmargin);

double median(uint64_t sorted_nums[], uint64_t len);

uint64_t rand_uint64();

#endif
