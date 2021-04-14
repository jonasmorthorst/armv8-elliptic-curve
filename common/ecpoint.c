#include <stdio.h>

#include "ecpoint.h"
#include "utils.h"

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem y, ef_elem z) {
  ec_point_lproj p;
	p.x = x;
  p.y = y;
	p.z = z;
	return p;
}

void ec_print_point_lproj_expr(ec_point_lproj p) {
  printf("x: ");
	ef_print_expr(p.x);
	printf(" y: ");
	ef_print_expr(p.y);
  printf(" z: ");
	ef_print_expr(p.z);
}

ec_point_lproj ec_rand_point_lproj() {
	return ec_create_point_lproj(ef_rand_elem(), ef_rand_elem(), ef_rand_elem());
}
