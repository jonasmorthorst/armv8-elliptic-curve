#include "extensionfield.h"

#ifndef ECPOINT_H
#define ECPOINT_H

#define EC_A00 0
#define EC_A01 0
#define EC_A10 1
#define EC_A11 0
#define EC_B00 0x2E6D944FA54DE7E5
#define EC_B01 0x59C8202CB9E6E0AE
#define EC_B10 0
#define EC_B11 0

typedef struct {
	ef_elem x, y, z;
} ec_point_lproj;

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem y, ef_elem z);

void ec_print_point_lproj_expr(ec_point_lproj p);

void ec_print_point_lproj_hex(ec_point_lproj p);

ec_point_lproj ec_rand_point_lproj();

ec_point_lproj ec_point_lproj_scalar_mull(uint64_t k, ec_point_lproj p);

ec_point_lproj ec_point_lproj_double_scalar_mull(uint64_t k, ec_point_lproj p1, uint64_t l, ec_point_lproj p2);

ec_point_lproj ec_point_lproj_neg(ec_point_lproj p);

ec_point_lproj ec_point_lproj_add(ec_point_lproj p1, ec_point_lproj p2);

ec_point_lproj ec_point_lproj_double(ec_point_lproj p);

#endif
