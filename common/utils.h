#include <arm_neon.h>

#include "basefield.h"
#include "extensionfield.h"
#include "ec.h"
#include "ec_scalarmull.h"

#ifndef UTILS_H
#define UTILS_H

void utils_init();

double average(uint64_t nums[], uint64_t len);

uint64_t compare_doubles(double a, double b, double errmargin);

uint64_t equal_poly64x2(poly64x2_t a, poly64x2_t b);

uint64_t equal_poly64x2x2(poly64x2x2_t a, poly64x2x2_t b);

uint64_t equal_ef_elem(ef_elem a, ef_elem b);

poly64x2x2_t concat_bf_poly(poly64x2_t p0, poly64x2_t p1);

double median(uint64_t sorted_nums[], uint64_t len);

uint64_t rand_uint64();

ec_point_laffine lproj_to_laffine(ec_point_lproj P);


#endif
