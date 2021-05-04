#include "ec_scalarmull.h"
#include <stdio.h>
#include "utils.h"
#include <stdlib.h>

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

ec_point_lproj ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64x1_t old_ptr, new_ptr, tmp, cond;
	cond[0]=1;
	int l = 64;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t c = 1-(decomp.k1[0]%2);
	decomp.k1[0] = decomp.k1[0]+c;

	// Compute recodings
	ec_naf naf_k1 = ec_to_naf(decomp.k1);
	ec_naf naf_k2 = ec_to_naf(decomp.k2);

	// Precomputation
	ec_point_laffine table[16];
	precompute(P, table);

	signed char k1_digit = naf_k1.val[l-1];
	uint64_t k1_sign = k1_digit < 0;

	signed char k2_digit = naf_k2.val[l-1];
	uint64_t k2_sign = k2_digit < 0;

	printf("k1 digit: %hhd | k2 digit: %hhd \n\n", k1_digit, k2_digit);

	ec_point_laffine P1 = table[abs(k1_digit)];
	ec_point_laffine P2 = ec_endo_laffine(table[abs(k2_digit)]);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_mixed(P1, ec_laffine_to_lproj(P2));

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(ec_double(Q)));

		k1_digit = naf_k1.val[i];
		k1_sign = k1_digit < 0;

		k2_digit = naf_k2.val[i];
		k2_sign = k2_digit < 0;

		printf("k1 digit: %hhd | k2 digit: %hhd \n\n", k1_digit, k2_digit);

		P1 = table[abs(k1_digit)];
		P2 = ec_endo_laffine(table[abs(k2_digit)]);

		P1_neg = ec_neg_laffine(P1);
		CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	uint64x2x2_t c_full = (uint64x2x2_t) {{{c, 0}, {0, 0}}};
	ec_point_lproj cP = ec_scalarmull_single(P, c_full);
	Q = ec_add(Q, ec_neg(cP));

	return Q;
}

void precompute(ec_point_laffine P, ec_point_laffine* table) {
	int end = 16;
	for (int i = 1; i < end; i++) {
		uint64x2x2_t k = (uint64x2x2_t) {{{i, 0}, {0, 0}}};
		table[i] = ec_lproj_to_laffine(ec_scalarmull_single(P, k));
	}
	// digitval[0]=1;
	// ec_point_laffine P_neg = ec_neg_laffine(P);
	// CMOV(tmp, decomp.k1_sign, digitval, P, P_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	//
	// uint64_t neg_Q_flag = decomp.k1_sign ^ decomp.k2_sign; //If sign is different they must be negated
	//
	// ec_point_laffine lookup[15]; //2^(w-1)-1 = 2^4 - 1 = 15
	// ec_point_lproj P2 = ec_double(ec_laffine_to_lproj(P));
	// ec_point_lproj P4 = ec_double(P2);
	// ec_point_lproj P8 = ec_double(P8);
	// lookup[0] = P;
	// lookup[1] = ec_lproj_to_laffine(P2);
	// lookup[2] = ec_lproj_to_laffine(ec_add_mixed(P, P2));
	// lookup[3] = ec_lproj_to_laffine(P4);
	// lookup[4] = ec_lproj_to_laffine(ec_add_mixed(P, P4));
	// lookup[5] = ec_lproj_to_laffine(ec_add_mixed(lookup[1], P4));
	// lookup[6] = ec_lproj_to_laffine(ec_add_mixed(lookup[2], P4));
	// lookup[7] = ec_lproj_to_laffine(P8);
	// lookup[8] = ec_lproj_to_laffine(ec_add_mixed(P, P8));
	// lookup[9] = ec_lproj_to_laffine(ec_add_mixed(lookup[1], P8));
	// lookup[10] = ec_lproj_to_laffine(ec_add_mixed(lookup[2], P8));
	// lookup[11] = ec_lproj_to_laffine(ec_add_mixed(lookup[3], P8));
	// lookup[12] = ec_lproj_to_laffine(ec_add_mixed(lookup[4], P8));
	// lookup[13] = ec_lproj_to_laffine(ec_add_mixed(lookup[5], P8));
	// lookup[14] = ec_lproj_to_laffine(ec_add_mixed(lookup[6], P8));
}

#define CEIL(A, B)                      (((A) - 1) / (B) + 1)
#define MASK(B)                         (((uint64_t)1 << (B)) - 1)

ec_naf ec_to_naf(uint64x2_t k) {
	uint64_t w = 5;
  //uint64_t m = 15; //2^(w-1);

	uint64_t order_len = 253;
	uint64_t l = CEIL(order_len, (w-1)); //64

	//0000 0001 1111
	uint64_t mask = MASK(w);

  ec_naf naf = { 0 };
  int i = 0;
  while (k[1] > 0 || i < l) {
		//int64_t k_i_temp = k[0]%(2*m)-m; // (k mod 16 only needs lower word)
    int64_t k_i_temp = (k[0] & mask) - ((int64_t)1 << (w - 1)); // (k mod 16 only needs lower word)

    // Extract 00001111
    //char digit = k_i_temp & 15;

    // If i is odd - store in higher bits
    //digit = digit << 4*(i&1);
    naf.val[i] = k_i_temp;

		// Subtraction
		k[0] -=k_i_temp;

    // Division by 8 => >> 3
    // We need to shift with carry

    // First shift lower bits
    k[0] = k[0] >> (w-1);

    // Then higher bits with 3 bits carry
    uint64_t carries = k[1] & 15;

    // Can we use both output input MOV R0,R0,ROR 4
    uint64_t shift_res;
    asm ("ROR %[res], %[input], #4;"
      : [res] "=r" (shift_res)
      : [input] "r" (carries)
      );

    k[0] = k[0] | shift_res;

    // First shift higher bits
    k[1] = k[1] >> (w-1);

    i++;
  }

  return naf;
}

void ec_print_naf(ec_naf k) {
  for(int i = 63; i >= 0; i--) {
    printf(" %hhd ", k.val[i]);
  }

  printf("\n");
}
