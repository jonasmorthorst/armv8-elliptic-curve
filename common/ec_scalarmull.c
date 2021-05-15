#include "ec_scalarmull.h"
#include <stdio.h>
#include "utils.h"
#include <stdlib.h>
#include <string.h>

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

void linear_pass(ec_point_laffine *P1, ec_point_laffine *P2, ec_point_laffine* table, uint64_t index1, uint64_t index2, uint64_t l) {
	uint64x1_t old_ptr, new_ptr, tmp, digitval;
	ec_point_laffine item, P1_tmp, P2_tmp;

	for (uint64_t i = 0; i < l; i++) {
		digitval[0]=i;
		item = table[i];
		CMOV(tmp, index1, digitval, P1_tmp, item, old_ptr, new_ptr, typeof(ec_point_laffine));
		CMOV(tmp, index2, digitval, P2_tmp, item, old_ptr, new_ptr, typeof(ec_point_laffine));
	}

	*P1 = P1_tmp;
	*P2 = ec_endo_laffine(P2_tmp);
}

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64x1_t old_ptr, new_ptr, tmp, cond;
	cond[0]=1;
	int l = 128;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char naf_k1[l];
	signed char naf_k2[l];

	ec_to_naf(decomp.k1, 3, naf_k1);
	ec_to_naf(decomp.k2, 3, naf_k2);

	// Precomputation
	ec_point_laffine table[2];
	precompute_w3(P, table);

	signed char k1_digit = naf_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = naf_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_mixed(P1, ec_laffine_to_lproj(P2));

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(Q);

		k1_digit = naf_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = naf_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 2);

		P1_neg = ec_neg_laffine(P1);
		CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed(P1, Q);
	CMOV(tmp, c1, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed(ec_endo_laffine(P2), Q);
	CMOV(tmp, c2, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64x1_t old_ptr, new_ptr, tmp, cond;
	cond[0]=1;
	int l = 86;

	//128
	//254

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char naf_k1[l];
	signed char naf_k2[l];

	ec_to_naf(decomp.k1, 4, naf_k1);
	ec_to_naf(decomp.k2, 4, naf_k2);

	// Precomputation
	ec_point_laffine table[4];
	precompute_w4(P, table);

	signed char k1_digit = naf_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = naf_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 4);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(Q));

		k1_digit = naf_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = naf_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 4);

		P1_neg = ec_neg_laffine(P1);
		CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CMOV(tmp, c1, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CMOV(tmp, c2, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	// printf("%s\n", "P IN");
	// ec_print_hex_laffine(P);

	uint64x1_t old_ptr, new_ptr, tmp, cond;
	cond[0]=1;
	int l = 65;
	// int l = 86;
	//l=52
	//86

	ec_split_scalar decomp = ec_scalar_decomp(k);

	// decomp.k1_sign = 1;
	uint64_t zero = 0;

	// printf("k1_0: %lu\n", decomp.k1[0]);
	// printf("k1_1: %lu\n", decomp.k1[1]);
	// printf("k2_0 %lu\n", decomp.k2[0]);
	// printf("k2_1 %lu\n", decomp.k2[1]);
	// printf("k1 sign %lu\n", decomp.k1_sign);
	// printf("k2 sign %lu\n", decomp.k2_sign);

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// printf("c1: %lu\n", c1);
	// printf("c2: %lu\n", c2);
	//
	// printf("k1_0: %lu\n", decomp.k1[0]);
	// printf("k1_1: %lu\n", decomp.k1[1]);
	// printf("k2_0 %lu\n", decomp.k2[0]);
	// printf("k2_1 %lu\n", decomp.k2[1]);

	// Compute recodings
	signed char naf_k1[l];
	signed char naf_k2[l];

	ec_to_naf(decomp.k1, 5, naf_k1);
	ec_to_naf(decomp.k2, 5, naf_k2);

	// ec_print_naf(naf_k1, l);
	// ec_print_naf(naf_k2, l);

	// Precomputation
	ec_point_laffine table[16];
	precompute(P, table);

	// print_table(table);

	signed char k1_digit = naf_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = naf_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;


	// printf("k1_digit_sign %lu\n", k1_digit_sign);
	// printf("k2_digit_sign %lu\n", k2_digit_sign);
	// printf("k1_sign %lu\n", k1_sign);
	// printf("k2_sign %lu\n", k2_sign);
	//
	// printf("k1 digit: %hhd | k2 digit: %hhd \n\n", k1_digit, k2_digit);

	ec_point_laffine P1;
	ec_point_laffine P2;

	linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 8);

	// ec_point_laffine P1 = table[k2_val/2];
	// ec_point_laffine P2 = ec_endo_laffine(table[k2_val/2]);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

	// printf("P1 negated: %lu\n", ec_equal_point_laffine(P1, P1_neg));
	// printf("P2 negated: %lu\n", ec_equal_point_laffine(P2, P2_neg));

	// ec_print_hex_laffine(P1);

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {

		// printf("Q On curve: %lu\n", ec_is_on_curve(Q));

		// printf("i: %d\n", i);
                                   
		Q = ec_double(ec_double(ec_double(Q)));

		k1_digit = naf_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = naf_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		// printf("k1 digit: %hhd | k2 digit: %hhd \n", k1_digit, k2_digit);
		// printf("k1 sign: %lu | k2 sign: %lu \n", k2_sign, k2_sign);
		// printf("k1 val: %hhd | k2 val: %hhd \n\n", k1_val, k2_digit);

		linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 8);

		// printf("P1 On curve: %lu\n", ec_is_on_curve(ec_laffine_to_lproj(P1)));
		// printf("P2 On curve: %lu\n", ec_is_on_curve(ec_laffine_to_lproj(P2)));

		// printf("%s\n", "P1 NEG");

		// ec_print_hex_laffine(P1_neg);

		//First XOR scalar signs

		//Negate p1 by k1_sign
		P1_neg = ec_neg_laffine(P1);
		CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

		// printf("%s\n", "P1");
		// ec_print_hex_laffine(P1);

		// printf("P1 negated: %lu\n", ec_equal_point_laffine(P1, P1_neg));
		// printf("P2 negated: %lu\n", ec_equal_point_laffine(P2, P2_neg));

		Q = ec_double_then_addtwo(P1, P2, Q);
		//
		// printf("Q after iteration i=%d\n", i);
		// ec_print_hex(Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CMOV(tmp, c1, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CMOV(tmp, c2, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	// printf("Q On curve: %lu\n", ec_is_on_curve(Q));

	// printf("Returning Q \n");
	// ec_print_hex(Q);

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64x1_t old_ptr, new_ptr, tmp, cond;
	cond[0]=1;
	int l = 52;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char naf_k1[l];
	signed char naf_k2[l];

	ec_to_naf(decomp.k1, 6, naf_k1);
	ec_to_naf(decomp.k2, 6, naf_k2);

	// Precomputation
	ec_point_laffine table[16];
	precompute_w6(P, table);

	signed char k1_digit = naf_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = naf_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 16);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(ec_double(ec_double(Q))));

		k1_digit = naf_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = naf_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		linear_pass(&P1, &P2, table, k1_val/2, k2_val/2, 16);

		//Negate p1 by k1_sign
		P1_neg = ec_neg_laffine(P1);
		CMOV(tmp, k1_sign, cond, P1, P1_neg, old_ptr, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CMOV(tmp, k2_sign, cond, P2, P2_neg, old_ptr, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CMOV(tmp, c1, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CMOV(tmp, c2, cond, Q, Q_add_neg, old_ptr, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

void precompute_first(ec_point_laffine P, ec_point_laffine* table) {
	ec_point_lproj P2 = ec_double(ec_laffine_to_lproj(P));
	table[0] = P;
	table[1] = ec_lproj_to_laffine(ec_add_mixed(table[1], P2));
	table[2] = ec_lproj_to_laffine(ec_add_mixed(table[3], P2));
	table[3] = ec_lproj_to_laffine(ec_add_mixed(table[5], P2));
	table[4] = ec_lproj_to_laffine(ec_add_mixed(table[7], P2));
	table[5] = ec_lproj_to_laffine(ec_add_mixed(table[9], P2));
	table[6] = ec_lproj_to_laffine(ec_add_mixed(table[11], P2));
	table[7] = ec_lproj_to_laffine(ec_add_mixed(table[13], P2));
}

void precompute_w3(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj Pl = ec_laffine_to_lproj(P);
	ec_point_lproj P3 = ec_double_then_add(P, Pl);

	ef_elem inv_inputs[1] = {P3.z};
	ef_elem inv_outputs[1];
	ef_sim_inv(inv_inputs, inv_outputs, 1);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_mull(P3.x, inv_outputs[0]), ef_mull(P3.l, inv_outputs[0])};
}

void precompute_w4(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P5 = ec_double_then_add(P, P2);
	ec_point_lproj P7 = ec_double_then_add(P, P3);
	
	ef_elem inv_inputs[3] = {P3.z, P5.z, P7.z};
	ef_elem inv_outputs[3];
	ef_sim_inv(inv_inputs, inv_outputs, 3);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_mull(P3.x, inv_outputs[0]), ef_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_mull(P5.x, inv_outputs[1]), ef_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_mull(P7.x, inv_outputs[2]), ef_mull(P7.l, inv_outputs[2])};
}

void precompute(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P4 = ec_double(P2);
	ec_point_lproj P5 = ec_add_mixed_unchecked(P, P4);
	ec_point_lproj P6 = ec_double(P3);
	ec_point_lproj P7 = ec_add_mixed_unchecked(P, P6);
	ec_point_lproj P9 = ec_double_then_add(P, P4);
	ec_point_lproj P11 = ec_double_then_add(P, P5);
	ec_point_lproj P13 = ec_double_then_add(P, P6);
	ec_point_lproj P15 = ec_double_then_add(P, P7);

	ef_elem inv_inputs[7] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z};
	ef_elem inv_outputs[7];
	ef_sim_inv(inv_inputs, inv_outputs, 7);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_mull(P3.x, inv_outputs[0]), ef_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_mull(P5.x, inv_outputs[1]), ef_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_mull(P7.x, inv_outputs[2]), ef_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_mull(P9.x, inv_outputs[3]), ef_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_mull(P11.x, inv_outputs[4]), ef_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_mull(P13.x, inv_outputs[5]), ef_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_mull(P15.x, inv_outputs[6]), ef_mull(P15.l, inv_outputs[6])};
}

void precompute_w6(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P4 = ec_double(P2);
	ec_point_lproj P5 = ec_add_mixed_unchecked(P, P4);
	ec_point_lproj P6 = ec_double(P3);
	ec_point_lproj P7 = ec_add_mixed_unchecked(P, P6);
	ec_point_lproj P8 = ec_double(P4);
	ec_point_lproj P9 = ec_add_mixed_unchecked(P, P8);
	ec_point_lproj P10 = ec_double(P5);
	ec_point_lproj P11 = ec_add_mixed_unchecked(P, P10);
	ec_point_lproj P12 = ec_double(P6);
	ec_point_lproj P13 = ec_add_mixed_unchecked(P, P12);
	ec_point_lproj P14 = ec_double(P7);
	ec_point_lproj P15 = ec_add_mixed_unchecked(P, P14);
	ec_point_lproj P17 = ec_double_then_add(P, P8);
	ec_point_lproj P19 = ec_double_then_add(P, P9);
	ec_point_lproj P21 = ec_double_then_add(P, P10);
	ec_point_lproj P23 = ec_double_then_add(P, P11);
	ec_point_lproj P25 = ec_double_then_add(P, P12);
	ec_point_lproj P27 = ec_double_then_add(P, P13);
	ec_point_lproj P29 = ec_double_then_add(P, P14);
	ec_point_lproj P31 = ec_double_then_add(P, P15);

	ef_elem inv_inputs[15] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z, P17.z, P19.z, P21.z, P23.z, P25.z, P27.z, P29.z, P31.z};
	ef_elem inv_outputs[15];
	ef_sim_inv(inv_inputs, inv_outputs, 15);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_mull(P3.x, inv_outputs[0]), ef_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_mull(P5.x, inv_outputs[1]), ef_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_mull(P7.x, inv_outputs[2]), ef_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_mull(P9.x, inv_outputs[3]), ef_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_mull(P11.x, inv_outputs[4]), ef_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_mull(P13.x, inv_outputs[5]), ef_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_mull(P15.x, inv_outputs[6]), ef_mull(P15.l, inv_outputs[6])};
	table[8] = (ec_point_laffine) {ef_mull(P17.x, inv_outputs[7]), ef_mull(P17.l, inv_outputs[7])};
	table[9] = (ec_point_laffine) {ef_mull(P19.x, inv_outputs[8]), ef_mull(P19.l, inv_outputs[8])};
	table[10] = (ec_point_laffine) {ef_mull(P21.x, inv_outputs[9]), ef_mull(P21.l, inv_outputs[9])};
	table[11] = (ec_point_laffine) {ef_mull(P23.x, inv_outputs[10]), ef_mull(P23.l, inv_outputs[10])};
	table[12] = (ec_point_laffine) {ef_mull(P25.x, inv_outputs[11]), ef_mull(P25.l, inv_outputs[11])};
	table[13] = (ec_point_laffine) {ef_mull(P27.x, inv_outputs[12]), ef_mull(P27.l, inv_outputs[12])};
	table[14] = (ec_point_laffine) {ef_mull(P29.x, inv_outputs[13]), ef_mull(P29.l, inv_outputs[13])};
	table[15] = (ec_point_laffine) {ef_mull(P31.x, inv_outputs[14]), ef_mull(P31.l, inv_outputs[14])};
}

void print_table(ec_point_laffine* table) {
	int end = 16;
	for (int i = 1; i < end; i++) {
		printf("\n %d) \n", i);
    ec_print_hex_laffine(table[i]);
  }

  printf("\n");
}

#define CEIL(A, B)                      (((A) - 1) / (B) + 1)
#define MASK(B)                         (((uint64_t)1 << (B)) - 1)

void ec_to_naf(uint64x2_t k, uint64_t w, signed char naf[]) {
  //uint64_t m = 15; //2^(w-1);

	uint64_t order_len = 253;
	uint64_t l = CEIL(order_len, (w-1)); //64

	//0000 0001 1111
	uint64_t mask = MASK(w);

  int i = 0;
  while (k[1] > 0 || i < l) {
		//int64_t k_i_temp = k[0]%(2*m)-m; // (k mod 16 only needs lower word)
    int64_t k_i_temp = (k[0] & mask) - ((int64_t)1 << (w - 1)); // (k mod 16 only needs lower word)

    // Extract 00001111
    //char digit = k_i_temp & 15;

    // If i is odd - store in higher bits
    //digit = digit << 4*(i&1);
    naf[i] = k_i_temp;

		// Subtraction
		k[0] -=k_i_temp;

    // Division by 8 => >> 3
    // We need to shift with carry

    // First shift lower bits
    k[0] = k[0] >> (w-1);

    // Then higher bits with 4 bits carry
    uint64_t carries = k[1] & mask;

    // Can we use both output input MOV R0,R0,ROR 4
    uint64_t shift_res;
    asm ("ROR %[res], %[input], %[w];"
      : [res] "=r" (shift_res)
      : [input] "r" (carries), [w] "r" (w-1)
      );

    k[0] = k[0] | shift_res;

    // First shift higher bits
    k[1] = k[1] >> (w-1);

    i++;
  }

	naf[i] = k[0] & mask;
}

void ec_print_naf(signed char *naf, uint64_t l) {
  for(int i = l-1; i >= 0; i--) {
    printf(" %hhd ", naf[i]);
  }

  printf("\n");
}

#define types_copy(Y,X) memcpy(Y, X, sizeof(uint64_t)*2);

void bn_rsh(uint64_t *a, int size, int bits) {
        int i;
        uint64_t r, carry, shift, mask, *c;

        a += size - 1;

        /* Prepare the bit mask. */
        shift = 64 - bits;
        carry = 0;
        mask = MASK(bits);
        for (i = size - 1; i >= 0; i--, a--) {
                /* Get the needed least significant bits. */
                r = (*a) & mask;
                /* Shift left the operand. */
                *a = ((*a) >> bits) | (carry << shift);
                /* Update the carry. */
                carry = r;
        }
}

void scmul_wreg(signed char *naf, int *len, elt k, int n, int w) {
        int i, l;
        elt t;
        int64_t t0, mask;
        int64_t u_i;
        int64_t efe=(uint64_t)(-1), zero = 0x0, one = 0x1;

        mask = MASK(w);
        l = CEIL(n, (w - 1));

        types_copy(t, k);

        i = 0;
        for (i = 0; t[1] != 0; i++, naf++) {
                u_i = (t[0] & mask) - ((uint64_t)1 << (w - 1));
                t[0] -= u_i;
                *naf = u_i;
                bn_rsh(t, 2, w - 1);
        }
        for (; i < l; i++, naf++) {
                u_i = (t[0] & mask) - ((uint64_t)1 << (w - 1));
                t[0] -= u_i;
                *naf = u_i;
                t[0] >>= (w - 1);
        }

        *naf = t[0] & mask;
        *len = l + 1;
}
