#include "ec_scalarmull.h"
#include <stdio.h>

#define pow2to63 9223372036854775808U

#define SUBACC_128(a0, a1, c0, c1)\
	asm volatile("SUBS %0, %0, %2;" \
		 "SBC %1, %1, %3;"\
		: "+r" (c0), "+r" (c1)\
		: "r" (a0), "r" (a1)\
		);

#define ADDACC_128(a0, a1, c0, c1)\
	asm volatile("ADDS %0, %0, %2;"\
		 "ADC %1, %1, %3;"\
		 : "+r" (c0), "+r" (c1)\
		 : "r" (a0), "r" (a1)\
		 );

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

ec_naf ec_to_naf(poly64x2_t k) {
  char m = 8;

  ec_naf naf = { 0 };

  // naf.val[0] = 209; //11010001
  // naf.val[1] = 91;  //01011011
  // naf.val[2] = 11;  //00001011

  int i = 0;


  while (k[1] > 0 || k[0] > m) {
    int naf_index = i/2; // Rounds up
    int64_t k_i_temp = k[0]%(2*m)-m; // (k mod 16 only needs lower word)

    // Extract 00001111
    char digit = k_i_temp & 15;

    // If i is odd - store in higher bits
    if (i & 1) {
      digit = digit << 4;
      naf.val[naf_index] = naf.val[naf_index] | digit;
    } else {
      naf.val[naf_index] = digit;
    }

		// Subtraction / addition
		uint64_t tmp1;
		uint64_t zero = 0;
		if (k_i_temp < 0) {
			tmp1 = -1*k_i_temp;
			ADDACC_128(tmp1, zero, k[0], k[1]);
		} else {
			SUBACC_128(k_i_temp, zero, k[0], k[1]);
		}

    // Division by 8 => >> 3
    // We need to shift with carry

    // First shift lower bits
    k[0] = k[0] >> 3;

    // Then higher bits with 3 bits carry
    uint64_t carries = k[1] & 7;

    // Can we use both output input MOV R0,R0,ROR 3
    uint64_t shift_res;
    asm ("ROR %[res], %[input], #3;"
      : [res] "=r" (shift_res)
      : [input] "r" (carries)
      );

    k[0] = k[0] | shift_res;

    // First shift higher bits
    k[1] = k[1] >> 3;

    i++;
  }

  return naf;
}

void ec_print_naf(ec_naf k) {

  // 1011 0001

  //d2 0000 1011
  for(int i = 15; i >= 0; i--) {
    //00001111
    unsigned char c1 = 15;
    //11110000
    unsigned char c2 = 240;
    //00001000
    unsigned char c1_sign_mask = 8;
    //00000111
    unsigned char c1_val_mask = 7;

    char d2 = ((char)k.val[i] & c2) >> 4;
    char d2_signed = d2 & c1_sign_mask;

    if (d2_signed) {
      printf(" %hhd ", (char)(d2 | c2));
    } else {
      printf(" %hhd ", d2);
    }

    char d1 = k.val[i] & c1;
    char d1_signed = d1 & c1_sign_mask;

    if (d1_signed) {
      printf(" %hhd ", (char)(d1 | c2));
    } else {
      printf(" %hhd ", d1);
    }
  }

  printf("\n");
}
