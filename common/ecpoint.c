#include <stdio.h>

#include "ecpoint.h"
#include "utils.h"

ec_point_lproj ec_create_point_lproj(ef_elem x, ef_elem y, ef_elem z) {
  ec_point_lproj p;
	p.x = x;
  p.l = y;
	p.z = z;
	return p;
}

void ec_print_point_lproj_expr(ec_point_lproj p) {
  printf("x: ");
	ef_print_expr(p.x);
	printf(" y: ");
	ef_print_expr(p.l);
  printf(" z: ");
	ef_print_expr(p.z);
  printf("\n");
}

void ec_print_point_lproj_hex(ec_point_lproj p) {
  printf("x: ");
	ef_print_hex(p.x);
	printf(" y: ");
	ef_print_hex(p.l);
  printf(" z: ");
	ef_print_hex(p.z);
  printf("\n");
}

ec_point_lproj ec_rand_point_lproj() {
	return ec_create_point_lproj(ef_rand_elem(), ef_rand_elem(), ef_rand_elem());
}

// Algorithm 3.27
ec_point_lproj ec_point_lproj_scalar_mull(uint64_t k, ec_point_lproj p) {
  // Init q to infinity
  // TODO: Make helper methods
  poly64x2_t p_one = {1, 0};
  poly64x2_t p_zero = {0, 0};

  ef_elem ef_one = ef_create_elem(p_one, p_zero);
  ef_elem ef_zero = ef_create_elem(p_zero, p_zero);

  ec_point_lproj q = ec_create_point_lproj(ef_one, ef_one, ef_zero);

  // To avoid adding and doubling with infinity
  int is_infty = 1;

  while (k) {
    if (!is_infty) {
      p = ec_point_lproj_double(p);
    }

    if (k & 1) {
      if (is_infty) {
        q = p;
        is_infty = 0;
      } else {
        q = ec_point_lproj_add(q, p);
      }
    }

    k = k/2;
  }

  return q;
}

// ec_point_lproj ec_point_lproj_double_scalar_mull(uint64_t k, ec_point_lproj p1, uint64_t l, ec_point_lproj p2) {
//
// }

ec_point_lproj ec_point_lproj_neg(ec_point_lproj p) {
  p.l = ef_add(p.l, p.z);

  return p;
}

// Assumption: q != +-p
ec_point_lproj ec_point_lproj_add(ec_point_lproj p, ec_point_lproj q) {
  printf("p: \n");
  ec_print_point_lproj_expr(p);

  printf("\nq: \n");
  ec_print_point_lproj_expr(q);

  ef_elem u = ef_add(ef_mull(p.y, q.z), ef_mull(q.y, p.z));

  ef_elem v = ef_square(ef_add(ef_mull(p.x, q.z), ef_mull(q.x, p.z)));

  ec_point_lproj r;
  r.x = ef_mull(ef_mull(ef_mull(u, ef_mull(p.x, q.z)), ef_mull(q.x, p.z)), u);
  r.l = ef_add(ef_square(ef_add(ef_mull(u, ef_mull(q.x, p.z)), v)), ef_mull(ef_mull(ef_mull(u, v), q.z), ef_add(p.l, p.z)));
  r.z = ef_mull(ef_mull(ef_mull(u, v), p.z), p.z);

  return r;
}

// Assumption: p != -p
ec_point_lproj ec_point_lproj_double(ec_point_lproj p) {
  // TODO: Make constant
  //poly64x2_t _a = {4611686018427387904, 4611686018427387904};
  //poly64x2_t _b = {4611686018427387904, 4611686018427387904};
  //ef_elem A = ef_create_elem(_a, _b);

  ef_elem w = ef_add(ef_add(ef_square(p.l), ef_mull(p.l, p.z)), ef_mull((ef_elem) A, ef_square(p.z)));

  ec_point_lproj r;
  r.x = ef_square(w);
  r.z = ef_mull(w, ef_square(r.z));
  r.l = ef_add(ef_add(ef_square(ef_mull(p.x, p.z)), r.x), ef_add(ef_mull(w, ef_mull(p.l, p.z)), r.z));

  return r;
}
