#include "ec_scalarmull.h"

// Algorithm 3.27
ec_point_lproj ec_scalarmull_single(ec_point_lproj P, poly64x2x2_t k) {
  // Init q to infinity
  // TODO: Make helper methods
  /*poly64x2_t p_one = {1, 0};
  poly64x2_t p_zero = {0, 0};

  ef_elem ef_one = ef_create_elem(p_one, p_zero);
  ef_elem ef_zero = ef_create_elem(p_zero, p_zero);

  ec_point_lproj Q = ec_create_point_lproj(ef_one, ef_one, ef_zero);

  // To avoid adding and doubling with infinity
  int is_infty = 1;

  while (k) {
    if (!is_infty) {
      P = ec_double(P);
    }

    if (k & 1) {
      if (is_infty) {
        Q = P;
        is_infty = 0;
      } else {
        Q = ec_add(Q, P);
      }
    }

    k = k/2;
  }

  return Q;*/
  return P;
}

ec_point_lproj ec_scalarmull_double(ec_point_lproj P1, poly64x2x2_t k1, ec_point_lproj P2, poly64x2x2_t k2) {
	return P1;
}
