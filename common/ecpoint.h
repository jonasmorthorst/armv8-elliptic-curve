#include "extensionfield.h"

#ifndef ECPOINT_H
#define ECPOINT_H

typedef struct {
	ef_elem x, y, z;
} ec_point_lproj;

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem y, ef_elem z);

void ec_print_point_lproj_expr(ec_point_lproj p);

void ec_print_point_lproj_hex(ec_point_lproj p);

ec_point_lproj ec_rand_point_lproj();

ec_point_lproj ec_point_lproj_scalar_mull(int k, ec_point_lproj p);

ec_point_lproj ec_point_lproj_double_scalar_mull(int k, ec_point_lproj p1, int l, ec_point_lproj p2);

ec_point_lproj ec_point_lproj_neg(ec_point_lproj p);

ec_point_lproj ec_point_lproj_add(ec_point_lproj p1, ec_point_lproj p2);

ec_point_lproj ec_point_lproj_double(ec_point_lproj p);

#endif
