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
  printf("\n");
}

void ec_print_point_lproj_hex(ec_point_lproj p) {
  printf("x: ");
	ef_print_hex(p.x);
	printf(" y: ");
	ef_print_hex(p.y);
  printf(" z: ");
	ef_print_hex(p.z);
  printf("\n");
}

ec_point_lproj ec_rand_point_lproj() {
	return ec_create_point_lproj(ef_rand_elem(), ef_rand_elem(), ef_rand_elem());
}

ec_point_lproj ec_point_lproj_neg(ec_point_lproj p) {
  p.y = ef_add(p.y, p.z);

  return p;
}

ec_point_lproj ec_point_lproj_add(ec_point_lproj p, ec_point_lproj q) {
  ef_elem u = ef_add(ef_mull(p.y, q.z), ef_mull(q.y, p.z));
  ef_elem v = ef_square(ef_add(ef_mull(p.x, q.z), ef_mull(q.x, p.z)));

  ec_point_lproj r;
  r.x = ef_mull(ef_mull(ef_mull(u, ef_mull(p.x, q.z)), ef_mull(q.x, p.z)), u);
  r.y = ef_add(ef_square(ef_add(ef_mull(u, ef_mull(q.x, p.z)), v)), ef_mull(ef_mull(ef_mull(u, v), q.z), ef_add(p.y, p.z)));
  r.z = ef_mull(ef_mull(ef_mull(u, v), p.z), p.z);

  return r;
}

ec_point_lproj ec_point_lproj_double(ec_point_lproj p) {
  // TODO: Make constant
  poly64x2_t _a = {4611686018427387904, 4611686018427387904};
  poly64x2_t _b = {4611686018427387904, 4611686018427387904};
  ef_elem A = ef_create_elem(_a, _b);

  ef_elem w = ef_add(ef_add(ef_square(p.y), ef_mull(p.y, p.z)), ef_mull(A, ef_square(p.z)));

  ec_point_lproj r;
  r.x = ef_square(w);
  r.z = ef_mull(w, ef_square(r.z));
  r.y = ef_add(ef_add(ef_square(ef_mull(p.x, p.z)), r.x), ef_add(ef_mull(w, ef_mull(p.y, p.z)), r.z));

  return r;
}
