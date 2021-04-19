#include "ec_scalarmull.h"
#include <stdio.h>

#define pow2to63 9223372036854775808U

// Algorithm 3.27
// k is a 253 bit number - therefore we use poly to represent it
ec_point_lproj ec_scalarmull_single(ec_point_lproj P, poly64x2x2_t k) {
  ec_point_lproj Q = (ec_point_lproj) INFTY;

  for(int i = 1; i >= 0; i--) {
    for(int j = 1; j >= 0; j--) {
      poly64_t c = pow2to63;

      for(int t = 63; t >= 0; t--) {
        Q = ec_double(Q);

        if (k.val[i][j] & c) {
          Q = ec_add(Q, P);
        }

        c = c/2;
      }
    }
  }

  return Q;
}

ec_point_lproj ec_scalarmull_double(ec_point_lproj P1, poly64x2x2_t k1, ec_point_lproj P2, poly64x2x2_t k2) {
	return P1;
}
