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

#define CMOV(tmpx1, cmpval, eqcondx1, old, new, old_ptrx1, new_ptrx1, point_type) \
	tmpx1[0] = cmpval & eqcondx1[0]; \
	tmpx1 = vceq_u64(tmpx1, eqcondx1); \
	old_ptrx1[0] = (uint64_t) &old; \
	new_ptrx1[0] = (uint64_t) &new; \
	tmpx1 = vbsl_u64(tmpx1, new_ptrx1, old_ptrx1); \
	old = *((point_type*) tmpx1[0]); \

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k) {
	ec_point_lproj new;
	uint64x1_t old_ptr, new_ptr, tmp, digitval;
	
	ec_split_scalar decomp = ec_scalar_decomp(k);
	ec_point_laffine Q = ec_endo_laffine(P);
	
	digitval[0]=1;
	ec_point_laffine P_neg = ec_neg_laffine(P);
	CMOV(tmp, decomp.k1_sign, digitval, P, P_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine Q_neg = ec_neg_laffine(Q);
	CMOV(tmp, decomp.k2_sign, digitval, Q, Q_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	
	ec_point_lproj R = (ec_point_lproj) INFTY;
	
	for(int i = 1; i >= 0; i--) {
		digitval[0] = pow2to63;
		for(int digit = 63; digit >= 0; digit--) {
			R = ec_double(R);

			new = ec_add_mixed(P, R);
			CMOV(tmp, decomp.k1[i], digitval, R, new, old_ptr, new_ptr, typeof(ec_point_lproj));
			
			new = ec_add_mixed(Q, R);
			CMOV(tmp, decomp.k2[i], digitval, R, new, old_ptr, new_ptr, typeof(ec_point_lproj));
			
			digitval[0] >>= 1;
		}
	}
	return ec_lproj_to_laffine(R);
}
/*
ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	//w=5 really only takes 4 bits at a time though
	ec_point_lproj new;
	uint64x1_t old_ptr, new_ptr, tmp,  ;
	
	ec_split_scalar decomp = ec_scalar_decomp(k);
	
	digitval[0]=1;
	ec_point_laffine P_neg = ec_neg_laffine(P);
	CMOV(tmp, decomp.k1_sign, digitval, P, P_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	
	uint64_t neg_Q_flag = decomp.k1_sign ^ decomp.k2_sign; //If sign is different they must be negated
	
	ec_point_laffine lookup[15]; //2^(w-1)-1 = 2^4 - 1 = 15
	ec_point_lproj P2 = ec_double(ec_laffine_to_lproj(P));
	ec_point_lproj P4 = ec_double(P2);
	ec_point_lproj P8 = ec_double(P8);
	lookup[0] = P;
	lookup[1] = ec_lproj_to_laffine(P2);
	lookup[2] = ec_lproj_to_laffine(ec_add_mixed(P, P2));
	lookup[3] = ec_lproj_to_laffine(P4);
	lookup[4] = ec_lproj_to_laffine(ec_add_mixed(P, P4));
	lookup[5] = ec_lproj_to_laffine(ec_add_mixed(lookup[1], P4));
	lookup[6] = ec_lproj_to_laffine(ec_add_mixed(lookup[2], P4));
	lookup[7] = ec_lproj_to_laffine(P8);
	lookup[8] = ec_lproj_to_laffine(ec_add_mixed(P, P8));
	lookup[9] = ec_lproj_to_laffine(ec_add_mixed(lookup[1], P8));
	lookup[10] = ec_lproj_to_laffine(ec_add_mixed(lookup[2], P8));
	lookup[11] = ec_lproj_to_laffine(ec_add_mixed(lookup[3], P8));
	lookup[12] = ec_lproj_to_laffine(ec_add_mixed(lookup[4], P8));
	lookup[13] = ec_lproj_to_laffine(ec_add_mixed(lookup[5], P8));
	lookup[14] = ec_lproj_to_laffine(ec_add_mixed(lookup[6], P8));
	
	ec_point_lproj R = (ec_point_lproj) INFTY;
	
	for(int i = 1; i >= 0; i--) {
		uint64_t windowmask = 0xF000000000000000;
		for(int j = 60; j >= 0; j -= 4) {
			R = ec_double(ec_double(ec_double(R)));
			uint64_t k1 = (decomp.k1[i] & windowmask) >> j;
			uint64_t k2 = (decomp.k2[i] & windowmask) >> j;
			
			
			windowmask >>= 4;
		}
	}
	return ec_lproj_to_laffine(R);
}*/
