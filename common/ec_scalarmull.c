#include "ec_scalarmull.h"
#include <stdio.h>

#define pow2to63 9223372036854775808U

// Algorithm 3.27
// k is a 253 bit number - therefore we use poly to represent it
ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k) {
	ec_point_lproj Q = (ec_point_lproj) INFTY;

	for(int i = 1; i >= 0; i--) {
		for(int j = 1; j >= 0; j--) {
			poly64_t c = pow2to63;

			for(int t = 63; t >= 0; t--) {
				Q = ec_double(Q);
				ec_point_lproj temp = ec_add_mixed(P, Q);

				uint64_t r0, val;
				asm ("ANDS %[val], %[k], %[c];"
					"CSEL %[r0], %[temp], %[q], NE;"
					: [val] "=r" ( val ), [r0] "=r" ( r0 )
					: [k] "r" (k.val[i][j]), [c] "r" (c), [q] "r" (&Q), [temp] "r" (&temp)
					);

				ec_point_lproj* ref = (ec_point_lproj*)r0;
				Q = *ref;

				c = c/2;
			}
		}
	}
	return Q;
}

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k) {
  ec_point_lproj Q = (ec_point_lproj) INFTY;

  for(int i = 1; i >= 0; i--) {
    for(int j = 1; j >= 0; j--) {
      poly64_t c = pow2to63;

      for(int t = 63; t >= 0; t--) {
        Q = ec_double(Q);

        if (k.val[i][j] & c) {
          Q = ec_add(P, Q);
        }

        c = c/2;
      }
    }
  }

  return Q;
}


#define pow2to64 {{{0, 1}, {0,0}}}

// Algorithm 3.48
ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l) {
  ec_point_lproj R = (ec_point_lproj) INFTY;

  for(int i = 1; i >= 0; i--) {
    for(int j = 1; j >= 0; j--) {
      R = ec_scalarmull_single_lproj(R, (uint64x2x2_t) pow2to64);

      ec_point_lproj k_i_P = ec_scalarmull_single_lproj(P, (uint64x2x2_t) {{{k.val[i][j], 0}, {0,0}}});
      ec_point_lproj l_i_Q = ec_scalarmull_single_lproj(Q, (uint64x2x2_t) {{{l.val[i][j], 0}, {0,0}}});

      R = ec_add(R, ec_add(k_i_P, l_i_Q));
    }
  }

  return R;
}

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k) {
	ec_split_scalar decomp = ec_scalar_decomp(k);
	ec_point_laffine Q = ec_endo_laffine(P);
	if(decomp.k1_sign) {
		P = ec_neg_laffine(P);
	}
	if(decomp.k2_sign) {
		Q = ec_neg_laffine(Q);
	}
	ec_point_lproj R = (ec_point_lproj) INFTY;
	ec_point_lproj new;
	uint64x1_t old_ptr, new_ptr, tmp, digitval;
	for(int i = 1; i >= 0; i--) {
		digitval[0] = pow2to63;
		for(int digit = 63; digit >= 0; digit--) {
			R = ec_double(R);

			new = ec_add_mixed(P, R);
			//cmov:
			tmp[0] = decomp.k1[i] & digitval[0];
			tmp = vceq_u64(tmp, digitval);
			old_ptr[0] = (uint64_t) &R;
			new_ptr[0] = (uint64_t) &new;
			tmp = vbsl_u64(tmp, new_ptr, old_ptr);
			R = *((ec_point_lproj*) tmp[0]);
			
			new = ec_add_mixed(Q, R);
			//cmov:
			tmp[0] = decomp.k2[i] & digitval[0];
			tmp = vceq_u64(tmp, digitval);
			old_ptr[0] = (uint64_t) &R;
			new_ptr[0] = (uint64_t) &new;
			tmp = vbsl_u64(tmp, new_ptr, old_ptr);
			R = *((ec_point_lproj*) tmp[0]);
			
			digitval[0] >>= 1;
		}
	}
	return ec_lproj_to_laffine(R);
}
