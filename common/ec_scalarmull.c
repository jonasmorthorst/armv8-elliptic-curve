#include "ec_scalarmull.h"
#include <stdlib.h>
#include <string.h>

#define pow2to63 9223372036854775808U
#define pow2to64 {{{0, 1}, {0,0}}}

// Algorithm 3.27
ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k) {
	ec_point_lproj Q = (ec_point_lproj) INFTY;

	for(int i = 1; i >= 0; i--) {
		for(int j = 1; j >= 0; j--) {
			poly64_t c = pow2to63;

			for(int t = 63; t >= 0; t--) {
				Q = ec_double(Q);
				ec_point_lproj temp = ec_add_mixed(P, Q);

				uint64_t r0, val;
				asm ("TST %[k], %[c];"
					"CSEL %[r0], %[temp], %[q], NE;"
					: [r0] "=r" ( r0 )
					: [k] "r" (k.val[i][j]), [c] "r" (c), [q] "r" (&Q), [temp] "r" (&temp)
					);
				ec_point_lproj* ref = (ec_point_lproj*)r0;
				Q = *ref;

				c >>= 1;
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

        c >>= 1;
      }
    }
  }

  return Q;
}

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

ec_point_laffine ec_scalarmull_single_endo_w3_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 65;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 3, rec_k1, l-1);
	reg_rec(decomp.k2, 3, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[2];
	precompute_w3(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));
	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_mixed(P1, ec_laffine_to_lproj(P2));

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(Q);

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));
		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));


	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w4_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 44;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 4, rec_k1, l-1);
	reg_rec(decomp.k2, 4, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[4];
	precompute_w4(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	lin_pass_w4(&P1, &P2, &table, k1_val/2, k2_val/2);
	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(Q));

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w4(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w5_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 33;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l], rec_k2[l];

	reg_rec(decomp.k1, 5, rec_k1, l-1);
	reg_rec(decomp.k2, 5, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[8];
	precompute_w5(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1, P2;
	lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(ec_double(Q)));

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w5(&P1, &P2, &table, k1_val/2, k2_val/2);

		P2 = ec_endo_laffine(P2);

		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	// Fix if c1 > 0
	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));


	// Fix if c2 > 0
	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

ec_point_laffine ec_scalarmull_single_endo_w6_randaccess(ec_point_laffine P, uint64x2x2_t k) {
	uint64_t new_ptr, con = 1;
	int l = 27;

	ec_split_scalar decomp = ec_scalar_decomp(k);

	uint64_t zero = 0;

	uint64_t c1 = 1-(decomp.k1[0]&1);
	decomp.k1[0] = decomp.k1[0]+c1;

	uint64_t c2 = 1-(decomp.k2[0]&1);
	decomp.k2[0] = decomp.k2[0]+c2;

	// Compute recodings
	signed char rec_k1[l];
	signed char rec_k2[l];

	reg_rec(decomp.k1, 6, rec_k1, l-1);
	reg_rec(decomp.k2, 6, rec_k2, l-1);

	// Precomputation
	ec_point_laffine table[16];
	precompute_w6(P, table);

	signed char k1_digit = rec_k1[l-1];
	uint64_t k1_digit_sign = ((unsigned char)k1_digit >> 7);
	signed char k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
	uint64_t k1_sign = k1_digit_sign^decomp.k1_sign;

	signed char k2_digit = rec_k2[l-1];
	uint64_t k2_digit_sign = (unsigned char)k2_digit >> 7;
	signed char k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
	uint64_t k2_sign = k2_digit_sign^decomp.k2_sign;

	ec_point_laffine P1;
	ec_point_laffine P2;

	lin_pass_w6(&P1, &P2, &table, k1_val/2, k2_val/2);
	P2 = ec_endo_laffine(P2);

	ec_point_laffine P1_neg = ec_neg_laffine(P1);
	CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

	ec_point_laffine P2_neg = ec_neg_laffine(P2);
	CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));


	ec_point_lproj Q = ec_add_laffine_unchecked(P1, P2);

	for(int i=l-2; i>=0; i--) {
		Q = ec_double(ec_double(ec_double(ec_double(Q))));

		k1_digit = rec_k1[i];
		k1_digit_sign = ((unsigned char)k1_digit >> 7);
		k1_val = (k1_digit^(zero - k1_digit_sign))+k1_digit_sign;
		k1_sign = k1_digit_sign^decomp.k1_sign;

		k2_digit = rec_k2[i];
		k2_digit_sign = ((unsigned char)k2_digit >> 7);
		k2_val = (k2_digit^(zero - k2_digit_sign))+k2_digit_sign;
		k2_sign = k2_digit_sign^decomp.k2_sign;

		lin_pass_w6(&P1, &P2, &table, k1_val/2, k2_val/2);
		P2 = ec_endo_laffine(P2);

		//Negate p1 by k1_sign
		P1_neg = ec_neg_laffine(P1);
		CSEL(k1_sign, con, P1, P1_neg, new_ptr, typeof(ec_point_laffine));

		P2_neg = ec_neg_laffine(P2);
		CSEL(k2_sign, con, P2, P2_neg, new_ptr, typeof(ec_point_laffine));

		Q = ec_double_then_addtwo(P1, P2, Q);
	}

	P1 = P;
	P1.l.val[0][0] ^= 1-decomp.k1_sign;
	ec_point_lproj Q_add_neg = ec_add_mixed_unchecked(P1, Q);
	CSEL(c1, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	P2 = P;
	P2.l.val[0][0] ^= 1-decomp.k2_sign;
	Q_add_neg = ec_add_mixed_unchecked(ec_endo_laffine(P2), Q);
	CSEL(c2, con, Q, Q_add_neg, new_ptr, typeof(ec_point_lproj));

	return ec_lproj_to_laffine(Q);
}

void precompute_w3(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj Pl = ec_laffine_to_lproj(P);
	ec_point_lproj P3 = ec_double_then_add(P, Pl);

	ef_intrl_elem P3Z_inv = ef_intrl_inv(P3.z);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, P3Z_inv), ef_intrl_mull(P3.l, P3Z_inv)};
}

void precompute_w4(ec_point_laffine P, ec_point_laffine table[]) {
	ec_point_lproj P2 = ec_double_mixed(P);
	ec_point_lproj P3 = ec_add_mixed_unchecked(P, P2);
	ec_point_lproj P5 = ec_double_then_add(P, P2);
	ec_point_lproj P7 = ec_double_then_add(P, P3);

	ef_intrl_elem inv_inputs[3] = {P3.z, P5.z, P7.z};
	ef_intrl_elem inv_outputs[3];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 3);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
}

void precompute_w5(ec_point_laffine P, ec_point_laffine table[]) {
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

	ef_intrl_elem inv_inputs[7] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z};
	ef_intrl_elem inv_outputs[7];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 7);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
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

	ef_intrl_elem inv_inputs[15] = {P3.z, P5.z, P7.z, P9.z, P11.z, P13.z, P15.z, P17.z, P19.z, P21.z, P23.z, P25.z, P27.z, P29.z, P31.z};
	ef_intrl_elem inv_outputs[15];
	ef_intrl_sim_inv(inv_inputs, inv_outputs, 15);

	table[0] = P;
	table[1] = (ec_point_laffine) {ef_intrl_mull(P3.x, inv_outputs[0]), ef_intrl_mull(P3.l, inv_outputs[0])};
	table[2] = (ec_point_laffine) {ef_intrl_mull(P5.x, inv_outputs[1]), ef_intrl_mull(P5.l, inv_outputs[1])};
	table[3] = (ec_point_laffine) {ef_intrl_mull(P7.x, inv_outputs[2]), ef_intrl_mull(P7.l, inv_outputs[2])};
	table[4] = (ec_point_laffine) {ef_intrl_mull(P9.x, inv_outputs[3]), ef_intrl_mull(P9.l, inv_outputs[3])};
	table[5] = (ec_point_laffine) {ef_intrl_mull(P11.x, inv_outputs[4]), ef_intrl_mull(P11.l, inv_outputs[4])};
	table[6] = (ec_point_laffine) {ef_intrl_mull(P13.x, inv_outputs[5]), ef_intrl_mull(P13.l, inv_outputs[5])};
	table[7] = (ec_point_laffine) {ef_intrl_mull(P15.x, inv_outputs[6]), ef_intrl_mull(P15.l, inv_outputs[6])};
	table[8] = (ec_point_laffine) {ef_intrl_mull(P17.x, inv_outputs[7]), ef_intrl_mull(P17.l, inv_outputs[7])};
	table[9] = (ec_point_laffine) {ef_intrl_mull(P19.x, inv_outputs[8]), ef_intrl_mull(P19.l, inv_outputs[8])};
	table[10] = (ec_point_laffine) {ef_intrl_mull(P21.x, inv_outputs[9]), ef_intrl_mull(P21.l, inv_outputs[9])};
	table[11] = (ec_point_laffine) {ef_intrl_mull(P23.x, inv_outputs[10]), ef_intrl_mull(P23.l, inv_outputs[10])};
	table[12] = (ec_point_laffine) {ef_intrl_mull(P25.x, inv_outputs[11]), ef_intrl_mull(P25.l, inv_outputs[11])};
	table[13] = (ec_point_laffine) {ef_intrl_mull(P27.x, inv_outputs[12]), ef_intrl_mull(P27.l, inv_outputs[12])};
	table[14] = (ec_point_laffine) {ef_intrl_mull(P29.x, inv_outputs[13]), ef_intrl_mull(P29.l, inv_outputs[13])};
	table[15] = (ec_point_laffine) {ef_intrl_mull(P31.x, inv_outputs[14]), ef_intrl_mull(P31.l, inv_outputs[14])};
}

#define MASK(B)                         (((uint64_t)1 << (B)) - 1)

void reg_rec(uint64x2_t k, uint64_t w, signed char rec[], uint64_t l) {
	uint64_t mask = MASK(w);

  int i = 0;
  while (k[1] > 0 || i < l) {
    int64_t k_i_temp = (k[0] & mask) - ((int64_t)1 << (w - 1)); // (k mod 16 only needs lower word)

    rec[i] = k_i_temp;

		k[0] -=k_i_temp;

    // First shift lower bits
    k[0] = k[0] >> (w-1);

    // Then higher bits with w-1 bits carry
    uint64_t carries = k[1] & mask;

    uint64_t shift_res;
    asm ("ROR %[res], %[input], %[s];"
      : [res] "=r" (shift_res)
      : [input] "r" (carries), [s] "r" (w-1)
      );

    k[0] = k[0] | shift_res;
    k[1] = k[1] >> (w-1);

    i++;
  }

	rec[i] = k[0] & mask;
}

void ec_print_rec(signed char *rec, uint64_t l) {
  for(int i = l-1; i >= 0; i--) {
    printf(" %hhd ", rec[i]);
  }

  printf("\n");
}
